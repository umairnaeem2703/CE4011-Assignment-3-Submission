# src/structural_validator.py

from banded_solver import UnstableStructureError
from parser import StructuralModel

class StructuralValidator:
    """
    Per-solve topology checks.

    This module identifies structural configurations that are fundamentally 
    unsolvable due to lack of boundary conditions or disconnected floating parts.
    """

    def __init__(self, model: StructuralModel):
        self.model = model

    def validate(self):
        """
        Run global topology checks. Raises UnstableStructureError if a 
        fatal condition is found. Prints warnings for non-fatal issues.
        """
        # ── Fatal Check 1: Zero Supports ────────────────────────────────
        # No supports at all guarantees rigid body motion.
        if not self.model.supports:
            raise UnstableStructureError(
                "No boundary conditions defined. "
                "The structure is entirely unsupported — "
                "rigid body modes exist. Global stiffness matrix K will be singular."
            )

        fatal_errors = []

        # ── Fatal Check 2: Global X-Restraint ───────────────────────────
        # Most 2D structures require at least one ux restraint to prevent sway.
        if not any(s.restrain_ux for s in self.model.supports.values()):
            fatal_errors.append(
                "No support restrains global x-translation. "
                "The structure can sway freely in the x-direction. "
                "Ensure at least one node is pinned or fixed."
            )

        # ── Connectivity Analysis ───────────────────────────────────────
        # Build adjacency list to find isolated sub-structures.
        adj = {}
        for el in self.model.elements.values():
            ni, nj = el.node_i.id, el.node_j.id
            adj.setdefault(ni, set()).add(nj)
            adj.setdefault(nj, set()).add(ni)

        # BFS to find connected components
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

        # ── Fatal Check 3: Floating Sub-structures ─────────────────────
        # Every connected part of the model must have at least one support.
        for component in components:
            has_support = any(n in self.model.supports for n in component)
            if not has_support:
                floating_members = sorted(
                    el.id for el in self.model.elements.values()
                    if el.node_i.id in component or el.node_j.id in component
                )
                fatal_errors.append(
                    f"Member(s) [{', '.join(floating_members)}] form a floating "
                    "sub-structure with no supports. Their DOFs will cause singular rows in K."
                )

        # ── Result Handling ─────────────────────────────────────────────
        if fatal_errors:
            raise UnstableStructureError(
                "\n".join(f"  ↳ {e}" for e in fatal_errors)
            )

        # ── Non-Fatal: Disconnected but supported ──────────────────────
        if len(components) > 1:
            print(f"INFO: {len(components)} disconnected sub-structures detected. "
                  "All parts are independently supported and solvable.")