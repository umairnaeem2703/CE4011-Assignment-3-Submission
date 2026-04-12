import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Dict, List

# ==========================================
# 1. DATA MODELS
# ==========================================

@dataclass
class Material:
    id: str
    E: float

@dataclass
class Section:
    id: str
    A: float
    I: float = 0.0

@dataclass
class Node:
    id: int
    x: float
    y: float
    dofs: List[int] = field(default_factory=list)

@dataclass
class Element:
    id: str
    type: str          
    node_i: Node       
    node_j: Node       
    material: Material
    section: Section
    release_start: bool = False
    release_end: bool = False

@dataclass
class Support:
    node: Node
    restrain_ux: bool = False
    restrain_uy: bool = False
    restrain_rz: bool = False

@dataclass
class PointLoad:
    node: Node
    fx: float = 0.0
    fy: float = 0.0
    mz: float = 0.0

@dataclass
class UDL:
    """Represents both legacy <udl> and new v2 <member_udl>"""
    element: Element
    wx: float = 0.0
    wy: float = 0.0
    is_local: bool = False
    fef_condition: str = "fixed-fixed"

@dataclass
class MemberPointLoad:
    element: Element
    position: float
    fx: float = 0.0
    fy: float = 0.0
    fef_condition: str = "fixed-fixed"

@dataclass
class LoadCase:
    id: str
    name: str = ""
    point_loads: List[PointLoad] = field(default_factory=list)
    udls: List[UDL] = field(default_factory=list)
    member_point_loads: List[MemberPointLoad] = field(default_factory=list)

@dataclass
class StructuralModel:
    name: str = "Untitled Model"
    materials: Dict[str, Material] = field(default_factory=dict)
    sections: Dict[str, Section] = field(default_factory=dict)
    nodes: Dict[int, Node] = field(default_factory=dict)
    elements: Dict[str, Element] = field(default_factory=dict)
    supports: Dict[int, Support] = field(default_factory=dict)
    load_cases: Dict[str, LoadCase] = field(default_factory=dict)

# ==========================================
# 2. THE PARSER LOGIC
# ==========================================

class XMLParser:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.model = StructuralModel()
        self.tree = ET.parse(filepath)
        self.root = self.tree.getroot()

    def parse(self) -> StructuralModel:
        """Executes the parsing pipeline and returns the populated model."""
        self.model.name = self.root.attrib.get('name', 'Untitled Model')
        
        self._parse_materials()
        self._parse_sections()
        self._parse_nodes()
        self._parse_elements()
        self._parse_boundaries()
        self._parse_loads()
        
        return self.model

    def _parse_materials(self):
        for mat in self.root.find('materials').findall('material'):
            m_id = mat.attrib['id']
            E = float(mat.attrib['E'])
            self.model.materials[m_id] = Material(id=m_id, E=E)

    def _parse_sections(self):
        for sec in self.root.find('sections').findall('section'):
            s_id = sec.attrib['id']
            A = float(sec.attrib.get('A', 0.0))
            I = float(sec.attrib.get('I', 0.0))
            self.model.sections[s_id] = Section(id=s_id, A=A, I=I)

    def _parse_nodes(self):
        for n in self.root.find('nodes').findall('node'):
            n_id = int(n.attrib['id'])
            x = float(n.attrib['x'])
            y = float(n.attrib['y'])
            self.model.nodes[n_id] = Node(id=n_id, x=x, y=y)

    def _parse_elements(self):
        elements_node = self.root.find('elements')
        
        for child in elements_node:
            if child.tag not in ['frame', 'truss']:
                continue
                
            e_id = child.attrib['id']
            node_i = self.model.nodes[int(child.attrib['node_i'])]
            node_j = self.model.nodes[int(child.attrib['node_j'])]
            material = self.model.materials[child.attrib['material']]
            section = self.model.sections[child.attrib['section']]
            
            release_start = False
            release_end = False
            
            releases = child.find('releases')
            if releases is not None:
                for rel in releases.findall('release'):
                    if rel.attrib.get('end') == 'i':
                        release_start = True
                    elif rel.attrib.get('end') == 'j':
                        release_end = True

            self.model.elements[e_id] = Element(
                id=e_id, type=child.tag,
                node_i=node_i, node_j=node_j,
                material=material, section=section,
                release_start=release_start, release_end=release_end
            )

    def _parse_boundaries(self):
        for sup in self.root.find('boundary_conditions').findall('support'):
            node = self.model.nodes[int(sup.attrib['node'])]
            sup_type = sup.attrib.get('type')
            
            ux = uy = rz = False
            
            if sup_type == 'fixed':
                ux = uy = rz = True
            elif sup_type == 'pin':
                ux = uy = True; rz = False
            elif sup_type == 'roller_x':
                ux = False; uy = True; rz = False
            elif sup_type == 'roller_y':
                ux = True; uy = False; rz = False
            else:
                ux = bool(int(sup.attrib.get('ux', 0)))
                uy = bool(int(sup.attrib.get('uy', 0)))
                rz = bool(int(sup.attrib.get('rz', 0)))

            self.model.supports[node.id] = Support(node, ux, uy, rz)

    def _parse_loads(self):
        loads_node = self.root.find('load_cases')
        if loads_node is None:
            return

        for lc_node in loads_node.findall('load_case'):
            lc_id = lc_node.attrib['id']
            lc_name = lc_node.attrib.get('name', '')
            lc = LoadCase(id=lc_id, name=lc_name)
            
            # Point loads
            for p_load in lc_node.findall('point_load'):
                node = self.model.nodes[int(p_load.attrib['node'])]
                fx = float(p_load.attrib.get('fx', 0.0))
                fy = float(p_load.attrib.get('fy', 0.0))
                mz = float(p_load.attrib.get('mz', 0.0))
                lc.point_loads.append(PointLoad(node, fx, fy, mz))
                
            # SCHEMA: member_udl
            for udl in lc_node.findall('member_udl'):
                element = self.model.elements[udl.attrib['element']]
                wx = float(udl.attrib.get('wx', 0.0))
                wy = float(udl.attrib.get('wy', 0.0))
                fef_cond = udl.attrib.get('fef_condition', 'fixed-fixed')
                lc.udls.append(UDL(element, wx, wy, False, fef_cond))

            # SCHEMA: member_point_load
            for mpl in lc_node.findall('member_point_load'):
                element = self.model.elements[mpl.attrib['element']]
                pos = float(mpl.attrib['position'])
                fx = float(mpl.attrib.get('fx', 0.0))
                fy = float(mpl.attrib.get('fy', 0.0))
                fef_cond = mpl.attrib.get('fef_condition', 'fixed-fixed')
                lc.member_point_loads.append(MemberPointLoad(element, pos, fx, fy, fef_cond))
                
            # Legacy: udl
            for udl in lc_node.findall('udl'):
                element = self.model.elements[udl.attrib['element']]
                wx = float(udl.attrib.get('wx', 0.0))
                wy = float(udl.attrib.get('wy', 0.0))
                is_local = udl.attrib.get('local', 'false').lower() == 'true'
                lc.udls.append(UDL(element, wx, wy, is_local, 'fixed-fixed'))
                
            self.model.load_cases[lc_id] = lc