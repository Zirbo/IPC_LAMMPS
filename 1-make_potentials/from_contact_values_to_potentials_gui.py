#! /usr/bin/python3

import sys, os
from PyQt5 import QtWidgets, uic

qtcreator_file  = "sources/from_contact_values_to_potentials_gui.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtcreator_file)


class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):

    SCRIPT = './sources/from_contact_values_to_potentials_gui.sh "{}" {} {}'

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        self.execute_button.clicked.connect(self.read_and_execute)


    def read_and_execute(self):
        try:
            # parse values
            model_name = self.model_box.toPlainText()
            model_type = self.ospc_box.currentText()
            enforce_ipc_geometry = self.ipc_box.isChecked()
            delta = float(self.delta_box.toPlainText())
            ecc_1 = float(self.ecc1_box.toPlainText())
            rad_1 = float(self.rad1_box.toPlainText())
            v_EE = float(self.vee_box.toPlainText())
            v_Ep1 = float(self.vep1_box.toPlainText())
            v_p1p1 = float(self.vp1p1_box.toPlainText())
            inputfile_content = "{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format(
                model_name, delta, ecc_1, rad_1, v_EE, v_Ep1, v_p1p1)

            if (model_type == "asymmetric"):
                ecc_2 = float(self.ecc2_box.toPlainText())
                rad_2 = float(self.rad2_box.toPlainText())
                v_Ep2 = float(self.vep2_box.toPlainText())
                v_p1p2 = float(self.vp1p2_box.toPlainText())
                v_p2p2 = float(self.vp2p2_box.toPlainText())
                inputfile_content += "{}\n{}\n{}\n{}\n{}\n".format(
                    ecc_2, rad_2, v_Ep2, v_p1p2, v_p2p2)

            # print recap text
            recap_string = "Model: {}\n".format(model_name)
            recap_string += "delta: {}\necc_1: {}, rad_1: {}\n".format(
                        delta, ecc_1, rad_1)
            if (model_type == "asymmetric"):
                recap_string += "ecc_2: {}, rad_2: {}\n".format(
                        ecc_2, rad_2)
            recap_string += "v_EE: {}, v_Ep1: {}, v_p1p1: {}\n".format(
                        v_EE, v_Ep1, v_p1p1)
            if (model_type == "asymmetric"):
                recap_string += "v_Ep2: {}, v_p1p2: {}, v_p2p2: {}\n".format(
                        v_Ep2, v_p1p2, v_p2p2)
            self.execute_window.setText(recap_string)

            # create inputfile
            with open("sources/inputfile", 'w') as inputfile:
                inputfile.write(inputfile_content)
            # run potential script
            print(str(self.SCRIPT.format(model_name, model_type,
                                         1 if enforce_ipc_geometry else 0)))
            res = os.popen(self.SCRIPT.format(model_name, model_type,
                                         1 if enforce_ipc_geometry else 0))
            print(res)
        except ValueError as e:
            self.show_msg_box("Input error", "Wrong value or empty input box.")
            self.execute_window.setText("Wrong value or empty input box.")




    def show_msg_box(self, title, text):
        err_msg = QtWidgets.QMessageBox()
        err_msg.setText(title)
        err_msg.setInformativeText(text)
        #err_msg.setStandardButtons(QtWidgets.QMessageBox.Retry)
        err_msg.setIcon(QtWidgets.QMessageBox.Warning)
        err_msg.show()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())
