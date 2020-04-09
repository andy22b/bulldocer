import numpy as np
import sys
from obspy.core import read
from matplotlib import dates as mdates
from matplotlib.backend_bases import key_press_handler
from PyQt5.QtCore import Qt

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(ApplicationWindow, self).__init__()
        self.p_pick_flag = False
        self.s_pick_flag = False
        st = read("IC.LSA.00.BHZ.M.2013.285.131629.SAC")
        tr = st[0]
        data = tr.data[:]
        times = tr.times()
        start_num = mdates.date2num(tr.stats.starttime.datetime)
        times_num = times/(3600. * 24) + start_num

        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)
        self.canvas = FigureCanvas(Figure(figsize=(5, 3)))
        self.canvas.setFocusPolicy(Qt.ClickFocus)
        self.canvas.setFocus()
        layout.addWidget(self.canvas)
        self.addToolBar(NavigationToolbar(self.canvas, self))

        self._static_ax = self.canvas.figure.subplots()
        self._static_ax.plot_date(times_num, data, "-")
        self._static_ax.

        self.canvas.mpl_connect('button_press_event', self.onclick)
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.canvas.mpl_connect('key_release_event', self.on_key_release)

    def on_key_press(self, event):
        print('you pressed', event.key)
        if event.key == 'p':
            self.p_pick_flag = True
        elif event.key == 's':
            self.s_pick_flag = True

    def on_key_release(self, event):
        print('you released', event.key)
        if event.key == 'p':
            self.p_pick_flag = False
        elif event.key == 's':
            self.s_pick_flag = False

    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))
        if self.p_pick_flag:
            print("Yay1")
        elif self.s_pick_flag:
            print("Yay2")


#def p_click(event):



if __name__ == "__main__":
    qapp = QtWidgets.QApplication(sys.argv)
    app = ApplicationWindow()
    app.show()
    qapp.exec_()