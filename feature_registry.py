# CCkit 功能说明：
"""Feature registry used by the CCkit desktop visualization shell.

The registry deliberately contains UI-neutral metadata.  A feature maps a
user-facing name to the script and preview workflow it owns; the Qt shell can
then provide a feature card without embedding every future modeling script in
one window class.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ViewerFeature:
    """Describe one visualization feature exposed by the desktop shell."""

    key: str
    name: str
    description: str
    script_path: Path | None
    available: bool = True
    card_title: str = "Parameters"


class FeatureRegistry:
    """Ordered collection of feature definitions keyed by stable identifiers."""

    def __init__(self, features: tuple[ViewerFeature, ...] = ()) -> None:
        self._features: dict[str, ViewerFeature] = {}
        for feature in features:
            self.register(feature)

    def register(self, feature: ViewerFeature) -> None:
        """Add or replace a feature while preserving insertion order."""

        if not feature.key or not feature.name:
            raise ValueError("viewer features require a key and display name")
        self._features[feature.key] = feature

    def get(self, key: str) -> ViewerFeature:
        """Return a feature or raise a clear error for an unknown key."""

        try:
            return self._features[key]
        except KeyError as exc:
            raise KeyError(f"unknown viewer feature: {key}") from exc

    def all(self) -> tuple[ViewerFeature, ...]:
        """Return features in their UI order."""

        return tuple(self._features.values())

    def available(self) -> tuple[ViewerFeature, ...]:
        """Return features whose implementation is currently connected."""

        return tuple(feature for feature in self.all() if feature.available)


def create_default_registry(
    repository_root: Path, *, core_directory: Path | None = None,
) -> FeatureRegistry:
    """Create the shell's built-in feature catalog.

    Bicrystal Builder, Surface Builder, Interface Builder and Wyckoff model are connected. Remaining entries are named
    extension points so future scripts can acquire a card without changing the
    viewer's top-level navigation or build-worker contract.
    """

    core_directory = core_directory or repository_root / "core"
    model_script = core_directory / "bicrystal" / "bicrystal_backend.py"
    return FeatureRegistry(
        (
            ViewerFeature(
                key="bicrystal",
                name="Bicrystal Builder",
                description=(
                    "Build bicrystal models and inspect symmetric-tilt grain boundaries."
                ),
                script_path=model_script,
                card_title="Bicrystal Builder parameters",
            ),
            ViewerFeature(
                key="wyckoff",
                name="Wyckoff model",
                description="Generate a crystal from a space group, lattice constants and Wyckoff positions.",
                script_path=core_directory / "wyckoff" / "wyckoff_model.py",
                card_title="Wyckoff model",
            ),
            ViewerFeature(
                key="surface",
                name="Surface Builder",
                description="Cut oriented surface slabs with prescribed in-plane directions, terminations and vacuum.",
                script_path=core_directory / "surface" / "surface_backend.py",
                card_title="Surface Builder",
            ),
            ViewerFeature(
                key="interface",
                name="Interface Builder",
                description="Build heterointerfaces with pymatgen or InterfaceMaster CSL/DSCL matching.",
                script_path=core_directory / "interface" / "interface_backend.py",
                available=True,
                card_title="Interface parameters",
            ),
            ViewerFeature(
                key="analysis",
                name="Structure Analysis",
                description="Reserved for future structure analysis workflows.",
                script_path=None,
                available=False,
                card_title="Analysis parameters",
            ),
        )
    )
