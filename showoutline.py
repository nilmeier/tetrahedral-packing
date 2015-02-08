# Show unit cell outline for a molecule.  Usage "runscript showoutline.py #0"
import UnitCell
d = UnitCell.show_unit_cell_dialog()
from chimera import specifier
m = specifier.evalSpec(arguments[0]).molecules()[0]
d.show_outline_model(m)
