# This Python file uses the following encoding: utf-8

#if__name__ == "__main__":
from PySide6.QtWidgets import QApplication
from dolomites import koil
from dolomites.koil_ui import *
import sys

# wind_a = koil.winding(24,2)
#
# wind_a.add_coil(koil.coil(3,22,1))
# wind_a.add_coil(koil.coil(4,21,1))
# wind_a.add_coil(koil.coil(15,10,1))
# wind_a.add_coil(koil.coil(16,9,1))
#
# print(wind_a.get_kw())
# print(wind_a.get_kw(2))
# print(wind_a.get_kw(10))
# print(wind_a.get_kw(14))
# print(wind_a.coils[0])

wind_b = koil.m_phase_winding()
# wind_b.compute_winding(3,12,5,single_layer=False)
wind_b.compute_winding(3,24,2,single_layer=True)
#wind_b.compute_winding(3,36,2,single_layer=True)
#wind_b.compute_winding(5,20,2,single_layer=True)

# #wind_b.compute_winding(3,9,4)
# #wind_b.compute_winding(3,3,1)
# print(wind_b.star)
# print(wind_b.windings[0].get_kw())
#
# print(wind_b.windings[0].coils)
# print(wind_b.windings[1].coils)
# print(wind_b.windings[2].coils)
#
# print(wind_b.windings[2].get_slot_matrix())
# print(wind_b.windings[2].get_slot_matrix('txt',True, name='a'))
# print(wind_b.windings[2].get_slot_matrix('txt',True, name='b'))
# print(wind_b.windings[2].get_slot_matrix('lua',True,name='a'))
# print(wind_b.windings[2].get_slot_matrix('getdp',True,name='a'))
# print(wind_b.windings[2].get_slot_matrix('getdp-2l',True,name='a'))
# print(wind_b.windings[0].get_slot_matrix('m-file',True,name='a'))
# print(wind_b.windings[1].get_slot_matrix('m-file',True,name='b'))
# print(wind_b.windings[2].get_slot_matrix('m-file',True,name='c'))
# print(wind_b.windings[0].get_getdp_circuit(id=100,name='a'))
# print(wind_b.windings[1].get_getdp_circuit(id=200,name='b'))
#
# print(wind_b.windings[2].get_getdp_circuit(id=300,name='c'))


# matplotlib.use('TkAgg')


# wind_b.calc_harmonics(1,cur=[1,-0.5,-0.50])
wind_b.calc_harmonics(1,cur=[0.866,0,-0.866])
x,y = wind_b.calc_mmf()

import matplotlib.pyplot as plt
plt.plot(x,y)
plt.show()


# exit()
if __name__ == "__main__":
    app = QApplication([])
    window = koil_widget()
    window.draw_winding(wind_b)
    window.show()
    sys.exit(app.exec())
