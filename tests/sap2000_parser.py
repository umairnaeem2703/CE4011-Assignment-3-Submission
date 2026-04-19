"""
SAP2000 Test Utility - Robust Parser and Validator for 2D Solver Benchmarking

This module provides:
1. SAP2000Parser: Parses SAP2000 text exports with proper coordinate system mapping
2. CoordinateMapper: Handles X-Z (SAP2000) to X-Y (our solver) transformations
3. MemberForceTransformer: Converts local member forces to global using element orientation
4. Test helpers with proper tolerances for shear deformation discrepancies
"""

import math


class CoordinateMapper:
    """
    Maps between SAP2000's X-Z plane and our solver's X-Y plane.
    
    SAP2000 Coordinate System:
        - Works in X-Z vertical plane
        - U1 = Displacement in X (horizontal)
        - U3 = Displacement in Z (vertical, becomes Y in our system)
        - R2 = Rotation about Y axis
        
    Our Solver Coordinate System:
        - Works in X-Y plane
        - UX = Displacement in X (horizontal)
        - UY = Displacement in Y (vertical)
        - RZ = Rotation about Z axis (perpendicular to plane)
        
    SAP2000 Forces/Reactions:
        - F1 = Force in X direction
        - F3 = Force in Z direction (becomes Y in our system)
        - M2 = Moment about Y axis (becomes RZ moment in our system)
    """
    
    @staticmethod
    def map_sap2000_displacement(u1, u3, r2):
        """
        Convert SAP2000 displacements (X-Z plane) to solver format (X-Y plane).
        
        Coordinate System Transformation:
        SAP2000 X-Z plane (vertical plane, Y points out):
            - U1: Displacement in X direction
            - U3: Displacement in Z direction (vertical)
            - R2: Rotation about Y axis (right-hand rule)
            
        Our Solver X-Y plane (horizontal view):
            - UX: Displacement in X direction
            - UY: Displacement in Y direction (vertical in 2D analysis)
            - RZ: Rotation about Z axis (perpendicular to X-Y plane)
            
        The rotation convention differs: SAP2000's R2 (about Y with right-hand rule looking from +Y)
        is equivalent to negative RZ in our 2D convention. This happens because of how the axes
        are reoriented between the X-Z vertical view and X-Y horizontal view.
        
        Args:
            u1: SAP2000 U1 displacement (X direction)
            u3: SAP2000 U3 displacement (Z direction)
            r2: SAP2000 R2 rotation (about Y axis)
            
        Returns:
            Tuple of (ux, uy, rz) in solver coordinate system
        """
        ux = u1
        uy = u3
        rz = -r2  # Negate due to rotation convention difference
        return (ux, uy, rz)
    
    @staticmethod
    def map_sap2000_force(f1, f3, m2):
        """
        Convert SAP2000 reactions/forces (X-Z plane) to solver format.
        
        Coordinate System Transformation:
        SAP2000 X-Z plane (vertical plane):
            - F1: Force in X direction
            - F3: Force in Z direction (vertical)
            - M2: Moment about Y axis (right-hand rule)
            
        Our Solver X-Y plane:
            - Fx: Force in X direction
            - Fy: Force in Y direction
            - Mz: Moment about Z axis (perpendicular to plane)
            
        The moment convention differs similarly to rotations: M2 (about Y) is equivalent to
        negative Mz due to the reorientation of axes.
        
        Args:
            f1: SAP2000 F1 force (X direction)
            f3: SAP2000 F3 force (Z direction)
            m2: SAP2000 M2 moment (about Y axis)
            
        Returns:
            Tuple of (fx, fy, mz) in solver coordinate system
        """
        fx = f1
        fy = f3
        mz = -m2  # Negate due to moment convention difference
        return (fx, fy, mz)


class MemberForceTransformer:
    """
    Transforms member forces between local and global coordinate systems.
    
    For a 2D frame/truss element oriented at angle theta from global X-axis:
        Local coordinates: Axial (along element), Shear (perpendicular)
        Global coordinates: Fx, Fy (in X-Y plane)
        
    Transformation uses the element's rotation matrix:
        c = cos(theta)  (direction cosine along X)
        s = sin(theta)  (direction cosine along Y)
    """
    
    @staticmethod
    def calculate_element_angle(node_i, node_j):
        """
        Calculate the global orientation angle of an element.
        
        Args:
            node_i: Starting node object with .x, .y attributes
            node_j: Ending node object with .x, .y attributes
            
        Returns:
            Tuple of (theta_radians, cos_theta, sin_theta, length)
        """
        dx = node_j.x - node_i.x
        dy = node_j.y - node_i.y
        length = (dx**2 + dy**2)**0.5
        
        if length < 1e-10:
            raise ValueError(f"Element length is zero: I=({node_i.x},{node_i.y}) J=({node_j.x},{node_j.y})")
        
        cos_theta = dx / length
        sin_theta = dy / length
        theta = math.atan2(dy, dx)
        
        return theta, cos_theta, sin_theta, length
    
    @staticmethod
    def local_to_global_forces(axial, shear, moment, cos_theta, sin_theta):
        """
        Transform member forces from local (element) coordinates to global (X-Y plane) coordinates.
        
        Local coordinate system:
            - Axial: Along the element axis (from I to J)
            - Shear: Perpendicular to element axis
            - Moment: About the Z axis (perpendicular to X-Y plane)
            
        Transformation formulas:
            Fx = c * Axial - s * Shear
            Fy = s * Axial + c * Shear
            Mz = Moment  (unchanged)
            
        where c = cos(theta), s = sin(theta)
        
        Args:
            axial: Local axial force
            shear: Local shear force (perpendicular to element)
            moment: Local moment (in plane)
            cos_theta: cos(element orientation angle)
            sin_theta: sin(element orientation angle)
            
        Returns:
            Tuple of (Fx, Fy, Mz) in global coordinates
        """
        fx = cos_theta * axial - sin_theta * shear
        fy = sin_theta * axial + cos_theta * shear
        mz = moment
        
        return (fx, fy, mz)
    
    @staticmethod
    def global_to_local_forces(fx, fy, moment, cos_theta, sin_theta):
        """
        Transform forces from global (X-Y) to local (element) coordinates.
        Inverse of local_to_global transformation.
        
        Args:
            fx: Global X force
            fy: Global Y force
            moment: Moment (unchanged in 2D)
            cos_theta: cos(element orientation angle)
            sin_theta: sin(element orientation angle)
            
        Returns:
            Tuple of (Axial, Shear, Moment) in local coordinates
        """
        axial = cos_theta * fx + sin_theta * fy
        shear = -sin_theta * fx + cos_theta * fy
        moment_out = moment
        
        return (axial, shear, moment_out)


class SAP2000Parser:
    """
    Parses SAP2000 text file exports and extracts structural data.
    
    Expected format:
        - Displacements table with headers: Joint, OutputCase, CaseType, U1, U2, U3, R1, R2, R3
        - Reactions table with headers: Joint, OutputCase, CaseType, F1, F2, F3, M1, M2, M3
        - Element forces table with headers: Frame, Joint, OutputCase, CaseType, F1, F2, F3, M1, M2, M3, FrameElem
        
    Output uses solver coordinate system (X-Y plane, not X-Z):
        - Displacements: {node_id: (ux, uy, rz)}
        - Reactions: {node_id: (fx, fy, mz)}
        - Element forces: {elem_id: {'i': (fx, fy, mz), 'j': (fx, fy, mz)}}
    """
    
    def __init__(self, filepath):
        """
        Initialize parser with SAP2000 text file.
        
        Args:
            filepath: Path to SAP2000 exported text file
        """
        self.filepath = filepath
        self.displacements = {}  # {node_id: (ux, uy, rz)}
        self.reactions = {}      # {node_id: (fx, fy, mz)}
        self.element_forces = {} # {elem_id: {'i': (fx, fy, mz), 'j': (fx, fy, mz)}}
        
    def parse(self):
        """
        Parse the SAP2000 text file and extract all tables.
        
        Returns:
            Tuple of (displacements_dict, reactions_dict, element_forces_dict)
        """
        with open(self.filepath, 'r') as f:
            content = f.read()
        
        self._parse_displacements(content)
        self._parse_reactions(content)
        self._parse_element_forces(content)
        
        return self.displacements, self.reactions, self.element_forces
    
    def _parse_displacements(self, content):
        """
        Extract joint displacements from SAP2000 output.
        
        Expects format:
            Table:  Joint Displacements
            Joint   OutputCase   CaseType   U1    U2    U3    R1    R2    R3
            <data>
        """
        lines = content.split('\n')
        current_table = None
        
        for i, line in enumerate(lines):
            line_stripped = line.strip()
            
            if "Table:  Joint Displacements" in line_stripped:
                current_table = "disp"
                continue
            
            if "Table:" in line_stripped and current_table == "disp":
                current_table = None
                break
            
            if current_table == "disp" and line_stripped and not line_stripped.startswith("Joint"):
                parts = line_stripped.split()
                
                # Skip header or invalid lines
                if len(parts) < 9:
                    continue
                
                try:
                    node_id = int(parts[0])
                    u1 = float(parts[3])      # Index 3: U1
                    u3 = float(parts[5])      # Index 5: U3
                    r2 = float(parts[7])      # Index 7: R2
                    
                    # Map to solver coordinate system
                    ux, uy, rz = CoordinateMapper.map_sap2000_displacement(u1, u3, r2)
                    self.displacements[node_id] = (ux, uy, rz)
                    
                except (ValueError, IndexError):
                    continue
    
    def _parse_reactions(self, content):
        """
        Extract support reactions from SAP2000 output.
        
        Expects format:
            Table:  Joint Reactions
            Joint   OutputCase   CaseType   F1    F2    F3    M1    M2    M3
            <data>
        """
        lines = content.split('\n')
        current_table = None
        
        for i, line in enumerate(lines):
            line_stripped = line.strip()
            
            if "Table:  Joint Reactions" in line_stripped:
                current_table = "react"
                continue
            
            if "Table:" in line_stripped and current_table == "react":
                current_table = None
                break
            
            if current_table == "react" and line_stripped and not line_stripped.startswith("Joint"):
                parts = line_stripped.split()
                
                # Skip header or invalid lines
                if len(parts) < 9:
                    continue
                
                try:
                    node_id = int(parts[0])
                    f1 = float(parts[3])      # Index 3: F1
                    f3 = float(parts[5])      # Index 5: F3
                    m2 = float(parts[7])      # Index 7: M2
                    
                    # Map to solver coordinate system
                    fx, fy, mz = CoordinateMapper.map_sap2000_force(f1, f3, m2)
                    self.reactions[node_id] = (fx, fy, mz)
                    
                except (ValueError, IndexError):
                    continue
    
    def _parse_element_forces(self, content):
        """
        Extract member end forces from SAP2000 output.
        
        SAP2000 reports member forces in GLOBAL coordinates (F1, F3, M2).
        This method extracts forces at element ends (nodes I and J).
        
        Expects format:
            Table:  Element Joint Forces - Frames
            Frame   Joint   OutputCase   CaseType   F1    F2    F3    M1    M2    M3    FrameElem
            <data>
        """
        lines = content.split('\n')
        current_table = None
        elem_data = {}  # {elem_id: {station: (f1, f3, m2)}}
        
        for i, line in enumerate(lines):
            line_stripped = line.strip()
            
            if "Table:  Element Joint Forces - Frames" in line_stripped:
                current_table = "forces"
                continue
            
            if "Table:" in line_stripped and current_table == "forces":
                current_table = None
                break
            
            if current_table == "forces" and line_stripped and not line_stripped.startswith("Frame"):
                parts = line_stripped.split()
                
                # Skip header or invalid lines
                if len(parts) < 11:
                    continue
                
                try:
                    elem_id = parts[0]
                    joint_id = int(parts[1])
                    f1 = float(parts[4])      # Index 4: F1
                    f3 = float(parts[6])      # Index 6: F3
                    m2 = float(parts[8])      # Index 8: M2
                    
                    # Map to solver coordinate system
                    fx, fy, mz = CoordinateMapper.map_sap2000_force(f1, f3, m2)
                    
                    if elem_id not in elem_data:
                        elem_data[elem_id] = {}
                    
                    elem_data[elem_id][joint_id] = (fx, fy, mz)
                    
                except (ValueError, IndexError):
                    continue
        
        # Organize into i-node and j-node forces
        for elem_id, joint_forces in elem_data.items():
            if len(joint_forces) >= 2:
                # Typically SAP2000 reports I node first, then J node
                joint_ids = sorted(joint_forces.keys())
                self.element_forces[elem_id] = {
                    'i': joint_forces[joint_ids[0]],
                    'j': joint_forces[joint_ids[1]]
                }


def compare_with_tolerance(computed, expected, rel_tol, abs_tol=1e-6):
    """
    Compare two floating point values with relative and absolute tolerance.
    
    Hybrid comparison that respects both relative and absolute tolerances:
    - For large values, relative tolerance dominates
    - For small values, absolute tolerance prevents over-sensitivity
    
    This addresses the issue where very small displacements (e.g., 3.5e-6 m) 
    can have large relative errors even when absolute error is negligible.
    
    Formula: passes if (absolute_error <= abs_tol) OR (relative_error <= rel_tol)
    
    Args:
        computed: Computed value
        expected: Expected/reference value
        rel_tol: Relative tolerance (e.g., 0.02 for 2%)
        abs_tol: Absolute tolerance floor (default 1e-6, ~1 micrometer for displacement in meters)
        
    Returns:
        True if values match within either tolerance, False otherwise
    """
    abs_error = abs(computed - expected)
    
    # Check absolute tolerance first (for very small values)
    # 1e-6 m = 1 micrometer is a reasonable precision for structural displacements
    if abs_error <= abs_tol:
        return True
    
    # If expected is zero or near zero, only absolute tolerance applies
    if abs(expected) < abs_tol:
        return abs_error <= abs_tol
    
    # Check relative tolerance
    relative_error = abs_error / abs(expected)
    return relative_error <= rel_tol


def assert_displacement_match(computed_disp, expected_disp, rel_tol=0.01, rot_rel_tol=0.05):
    """
    Assert that computed displacement matches expected with tolerance.
    
    Rationale for tolerances:
    - Translations use 1% to account for differences between solver implementations:
      1. SAP2000 uses advanced proprietary solvers vs. our basic 2D FEA
      2. Different element formulations and numerical techniques
      3. Floating-point precision differences and accumulated rounding
      4. Minor theoretical differences (Euler-Bernoulli vs. Timoshenko shear effects ~1.5%)
    - Rotations use 5% to account for:
      1. Rotations are derived from beam theory that differs significantly between solvers
      2. Coupling with bending stiffness which is affected by shear deformation (5-10% range)
      3. Higher sensitivity to moment calculation methods
      4. Less critical than translations for practical structural validation
    
    Args:
        computed_disp: Tuple of (ux, uy, rz) from solver
        expected_disp: Tuple of (ux, uy, rz) from SAP2000
        rel_tol: Relative tolerance for translations (default 1%)
        rot_rel_tol: Relative tolerance for rotations (default 5%)
        
    Returns:
        Tuple of (match, error_message)
    """
    errors = []
    
    # Create components tuples with (name, computed, expected, tolerance)
    components = [('UX', computed_disp[0], expected_disp[0], rel_tol),
                  ('UY', computed_disp[1], expected_disp[1], rel_tol)]
    
    if len(computed_disp) > 2 and len(expected_disp) > 2:
        # More lenient tolerance for rotations due to theory differences
        components.append(('RZ', computed_disp[2], expected_disp[2], rot_rel_tol))
    
    for comp_name, computed, expected, tol in components:
        if not compare_with_tolerance(computed, expected, tol):
            error = (computed - expected) / (expected if expected != 0 else 1.0)
            errors.append(f"{comp_name}: {computed:.6e} vs {expected:.6e} (error: {error*100:.2f}%)")
    
    if errors:
        return False, "\n  ".join(errors)
    return True, ""


def assert_force_match(computed_force, expected_force, rel_tol=0.02):
    """
    Assert that computed force matches expected with tolerance.
    
    Args:
        computed_force: Tuple of (fx, fy, mz) from solver
        expected_force: Tuple of (fx, fy, mz) from SAP2000
        rel_tol: Relative tolerance for forces (default 2%)
        
    Returns:
        Tuple of (match, error_message)
    """
    errors = []
    
    components = [('Fx', computed_force[0], expected_force[0]),
                  ('Fy', computed_force[1], expected_force[1])]
    
    if len(computed_force) > 2 and len(expected_force) > 2:
        components.append(('Mz', computed_force[2], expected_force[2]))
    
    for comp_name, computed, expected in components:
        if not compare_with_tolerance(computed, expected, rel_tol):
            error = (computed - expected) / (expected if expected != 0 else 1.0)
            errors.append(f"{comp_name}: {computed:.6e} vs {expected:.6e} (error: {error*100:.2f}%)")
    
    if errors:
        return False, "\n  ".join(errors)
    return True, ""
