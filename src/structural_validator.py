# src/structural_validator.py

from banded_solver import UnstableStructureError
from parser import StructuralModel

class StructuralValidator:
    """Validates structural topology before solving (boundary conditions, connectivity)."""

    def __init__(self, model: StructuralModel):
        self.model = model

    def validate(self):
        """Checks for fatal structural issues: missing supports, unrestrained translations, floating parts."""
        # Check for zero supports
        if not self.model.supports:
            raise UnstableStructureError(
                "No boundary conditions defined. Structure is entirely unsupported."
            )

        fatal_errors = []

        # Check for global x-restraint
        if not any(s.restrain_ux for s in self.model.supports.values()):
            fatal_errors.append(
                "No support restrains global x-translation. "
                "Ensure at least one node is pinned or fixed."
            )

        # Build adjacency graph and find connected components
        adj = {}
        for el in self.model.elements.values():
            ni, nj = el.node_i.id, el.node_j.id
            adj.setdefault(ni, set()).add(nj)
            adj.setdefault(nj, set()).add(ni)

        unvisited = set(adj)
        components = []
        while unvisited:
            start = next(iter(unvisited))
            queue, component = [start], set()
            while queue:
                node = queue.pop(0)
                if node in component:
                    continue
                component.add(node)
                queue.extend(n for n in adj.get(node, []) if n not in component)
            unvisited -= component
            components.append(frozenset(component))

        # Check for floating sub-structures
        for component in components:
            has_support = any(n in self.model.supports for n in component)
            if not has_support:
                floating_members = sorted(
                    el.id for el in self.model.elements.values()
                    if el.node_i.id in component or el.node_j.id in component
                )
                fatal_errors.append(
                    f"Floating elements [{', '.join(floating_members)}] have no supports."
                )

        if fatal_errors:
            raise UnstableStructureError("\n".join(f"  ↳ {e}" for e in fatal_errors))

        if len(components) > 1:
            print(f"INFO: {len(components)} disconnected sub-structures detected (all supported).")