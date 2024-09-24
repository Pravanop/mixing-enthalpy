from pymatgen.core import Composition
from pymatgen.entries import Entry


class PDEntryLocal(Entry):
    """
    A local version of PDEntry that allows for additional attributes to be stored.
    """
    def __init__(
        self,
        composition: Composition,
        energy: float,
        name: str = None,
        attribute: object = None,
    ):
        super().__init__(composition, energy)
        self.name = name
        self.attribute = attribute

    def __repr__(self):
        name = ""
        if self.name:
            name = f" ({self.name})"
        return f"{type(self).__name__} : {self.composition}{name} with energy = {self.energy:.4f}"

    def update_energy(self, energy: float):
        self._energy = energy

    @property
    def energy(self) -> float:
        """The entry's energy."""
        return self._energy

    def as_dict(self):
        """Get MSONable dict representation of PDEntry."""
        return super().as_dict() | {"name": self.name, "attribute": self.attribute}

    @classmethod
    def from_dict(cls, dct: dict):
        """
        Args:
            dct (dict): dictionary representation of PDEntry.

        Returns:
            PDEntry
        """
        return cls(
            composition=Composition(dct["composition"]),
            energy=dct["energy"],
            name=dct.get("name"),
            attribute=dct.get("attribute"),
        )
