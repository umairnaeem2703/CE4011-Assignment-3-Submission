import xml.etree.ElementTree as ET
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List

# ==========================================
# 1. DATA MODELS
# ==========================================

@dataclass
class Material:
    id: str
    E: float
    alpha: float = 0.0  # Coefficient of thermal expansion

@dataclass
class Section:
    id: str
    A: float
    I: float = 0.0
    d: float = 0.0      # Section depth

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

# --- LOAD OOP HIERARCHY (From Diagram) ---

class Load(ABC):
    @abstractmethod
    def NodalLoads(self) -> list:
        pass

    @abstractmethod
    def FEF(self, fef_condition: str, L: float) -> list:
        pass

@dataclass
class NodalLoad(Load):
    node: Node
    fx: float = 0.0
    fy: float = 0.0
    mz: float = 0.0

    def NodalLoads(self) -> list:
        return [self.fx, self.fy, self.mz]

    def FEF(self, fef_condition: str, L: float) -> list:
        return [[0.0] for _ in range(6)]

@dataclass
class MemberLoad(Load, ABC):
    element: Element
    
    def NodalLoads(self) -> list:
        return [0.0, 0.0, 0.0]

@dataclass
class PointLoad(MemberLoad):
    position: float
    fx: float = 0.0
    fy: float = 0.0

    def FEF(self, fef_condition: str, L: float) -> list:
        fef = [[0.0] for _ in range(6)]
        a = self.position
        b = L - a
        # CRITICAL FIX: Use magnitude of load in FEF formulas (textbook convention)
        # The sign indicates load direction; formulas compute reaction magnitudes
        P = abs(self.fy)  
        fx_load = self.fx

        # Axial FEF: preserve sign for compression/tension
        fef[0][0] = fx_load * b / L
        fef[3][0] = fx_load * a / L

        # Transverse FEF: use magnitude, formulas give correct reaction magnitude
        if fef_condition == "fixed-fixed":
            fef[1][0] = P * (b**2) * (3*a + b) / (L**3)
            fef[2][0] = P * a * (b**2) / (L**2)
            fef[4][0] = P * (a**2) * (3*b + a) / (L**3)
            fef[5][0] = -P * (a**2) * b / (L**2)
        elif fef_condition == "pin-fixed":
            Mz_j = P * a * b * (L + a) / (2 * L**2)
            fef[1][0] = (P * b - Mz_j) / L
            fef[2][0] = 0.0
            fef[4][0] = (P * a + Mz_j) / L
            fef[5][0] = Mz_j
        elif fef_condition == "fixed-pin":
            Mz_i = P * a * b * (L + b) / (2 * L**2)
            fef[1][0] = (P * b + Mz_i) / L
            fef[2][0] = Mz_i
            fef[4][0] = (P * a - Mz_i) / L
            fef[5][0] = 0.0
        elif fef_condition == "pin-pin":
            fef[1][0] = P * b / L
            fef[2][0] = 0.0
            fef[4][0] = P * a / L
            fef[5][0] = 0.0
            
        return fef

@dataclass
class UniformlyDL(MemberLoad):
    wx: float = 0.0
    wy: float = 0.0
    is_local: bool = False

    def FEF(self, fef_condition: str, L: float) -> list:
        fef = [[0.0] for _ in range(6)]
        # CRITICAL FIX: Use magnitude of distributed load in FEF formulas
        # Sign indicates direction; formulas compute reaction magnitudes
        w = abs(self.wy)
        wx = self.wx  # Preserve sign for axial consistency

        # Axial FEF
        fef[0][0] = wx * L / 2.0
        fef[3][0] = wx * L / 2.0

        # Transverse FEF: use magnitude, formulas give correct reaction magnitude
        if fef_condition == "fixed-fixed":
            fef[1][0] = w * L / 2.0
            fef[2][0] = w * (L**2) / 12.0
            fef[4][0] = w * L / 2.0
            fef[5][0] = -w * (L**2) / 12.0
        elif fef_condition == "pin-fixed":
            fef[1][0] = (3.0 / 8.0) * w * L
            fef[2][0] = 0.0
            fef[4][0] = (5.0 / 8.0) * w * L
            fef[5][0] = w * (L**2) / 8.0
        elif fef_condition == "fixed-pin":
            fef[1][0] = (5.0 / 8.0) * w * L
            fef[2][0] = w * (L**2) / 8.0
            fef[4][0] = (3.0 / 8.0) * w * L
            fef[5][0] = 0.0
        elif fef_condition == "pin-pin":
            fef[1][0] = w * L / 2.0
            fef[2][0] = 0.0
            fef[4][0] = w * L / 2.0
            fef[5][0] = 0.0
            
        return fef

@dataclass
class TemperatureL(MemberLoad):
    Tu: float = 0.0
    Tb: float = 0.0

    def FEF(self, fef_condition: str, L: float) -> list:
        fef = [[0.0] for _ in range(6)]
        
        E = self.element.material.E
        alpha = self.element.material.alpha
        A = self.element.section.A
        I = self.element.section.I
        d = self.element.section.d
        
        # Decompose trapezoidal profile
        T_uniform = (self.Tu + self.Tb) / 2.0
        T_grad = self.Tu - self.Tb
        
        # Axial FEF (independent of bending releases)
        fef[0][0] = -alpha * T_uniform * E * A
        fef[3][0] =  alpha * T_uniform * E * A
        
        if self.element.type == 'truss':
            return fef
            
        # Fixed-Fixed Base Thermal Moments
        M_zi = -(alpha * T_grad / d) * E * I if d != 0 else 0.0
        M_zj =  (alpha * T_grad / d) * E * I if d != 0 else 0.0
        
        # Adjust moments and induced shears based on end releases
        if fef_condition == "fixed-fixed":
            fef[2][0] = M_zi
            fef[5][0] = M_zj
        elif fef_condition == "pin-fixed":
            fef[2][0] = 0.0
            fef[5][0] = M_zj - 0.5 * M_zi
            fef[1][0] = -1.5 * M_zi / L
            fef[4][0] =  1.5 * M_zi / L
        elif fef_condition == "fixed-pin":
            fef[2][0] = M_zi - 0.5 * M_zj
            fef[5][0] = 0.0
            fef[1][0] =  1.5 * M_zj / L
            fef[4][0] = -1.5 * M_zj / L
        elif fef_condition == "pin-pin":
            fef[2][0] = 0.0
            fef[5][0] = 0.0
            
        return fef

@dataclass
class LoadCase:
    id: str
    name: str = ""
    loads: List[Load] = field(default_factory=list)

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
            alpha = float(mat.attrib.get('alpha', 0.0))
            self.model.materials[m_id] = Material(id=m_id, E=E, alpha=alpha)

    def _parse_sections(self):
        for sec in self.root.find('sections').findall('section'):
            s_id = sec.attrib['id']
            A = float(sec.attrib.get('A', 0.0))
            I = float(sec.attrib.get('I', 0.0))
            d = float(sec.attrib.get('d', 0.0))
            self.model.sections[s_id] = Section(id=s_id, A=A, I=I, d=d)

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
            
            # Point loads -> NodalLoad
            for p_load in lc_node.findall('point_load'):
                node = self.model.nodes[int(p_load.attrib['node'])]
                fx = float(p_load.attrib.get('fx', 0.0))
                fy = float(p_load.attrib.get('fy', 0.0))
                mz = float(p_load.attrib.get('mz', 0.0))
                lc.loads.append(NodalLoad(node, fx, fy, mz))
                
            # SCHEMA: member_udl -> UniformlyDL
            for udl in lc_node.findall('member_udl'):
                element = self.model.elements[udl.attrib['element']]
                wx = float(udl.attrib.get('wx', 0.0))
                wy = float(udl.attrib.get('wy', 0.0))
                lc.loads.append(UniformlyDL(element, wx, wy, False))

            # SCHEMA: member_point_load -> PointLoad
            for mpl in lc_node.findall('member_point_load'):
                element = self.model.elements[mpl.attrib['element']]
                pos = float(mpl.attrib['position'])
                fx = float(mpl.attrib.get('fx', 0.0))
                fy = float(mpl.attrib.get('fy', 0.0))
                lc.loads.append(PointLoad(element, pos, fx, fy))
                
            # SCHEMA: temperature_load -> TemperatureL
            for tload in lc_node.findall('temperature_load'):
                element = self.model.elements[tload.attrib['element']]
                load_type = tload.attrib.get('type', 'uniform')
                
                if load_type == 'uniform':
                    dT = float(tload.attrib.get('delta_T', 0.0))
                    lc.loads.append(TemperatureL(element, Tu=dT, Tb=dT))
                elif load_type == 'gradient':
                    if element.type == 'truss':
                        raise ValueError(f"Truss element {element.id} cannot take gradient temperature loads.")
                    dT = float(tload.attrib.get('delta_T', 0.0))
                    lc.loads.append(TemperatureL(element, Tu=dT/2.0, Tb=-dT/2.0))
                elif load_type == 'combined':
                    if element.type == 'truss':
                        raise ValueError(f"Truss element {element.id} cannot take combined temperature loads.")
                    ttop = float(tload.attrib.get('T_top', 0.0))
                    tbot = float(tload.attrib.get('T_bottom', 0.0))
                    lc.loads.append(TemperatureL(element, Tu=ttop, Tb=tbot))

            # Legacy: udl -> UniformlyDL
            for udl in lc_node.findall('udl'):
                element = self.model.elements[udl.attrib['element']]
                wx = float(udl.attrib.get('wx', 0.0))
                wy = float(udl.attrib.get('wy', 0.0))
                is_local = udl.attrib.get('local', 'false').lower() == 'true'
                lc.loads.append(UniformlyDL(element, wx, wy, is_local))
                
            self.model.load_cases[lc_id] = lc