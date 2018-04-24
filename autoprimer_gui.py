#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 10:19:18 2018

@author: HCC604
"""
import sys

from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QLabel, QPushButton, QAction, QLineEdit, QMessageBox, QComboBox, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot, QThread

import autoprimer

class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'AutoPrimer GUI'
        self.left = 200
        self.top = 200
        self.width = 300
        self.height = 200
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # gene name input box
        self.label1 = QLabel(self)
        self.label1.move(20, 20)
        self.label1.setText('Gene Name:')
        self.textbox = QLineEdit(self)
        self.textbox.move(80, 20)
        self.textbox.resize(80, 30)


        # transcript combo box
        self.combobox = QComboBox(self)
        self.combobox.move(180, 20)

        # server combo box
        self.label2 = QLabel(self)
        self.label2.move(20, 80)
        self.label2.setText('SNPfree Server:')
        self.combobox_server = QComboBox(self)
        self.combobox_server.move(180, 80)
        server_status = autoprimer.check_servers()
        if str(server_status['hku_status']) == '200':
            self.combobox_server.addItem('HKU')
        if str(server_status['pyn_status']) == '200':
            self.combobox_server.addItem('PYN')
        if str(server_status['cyk_status']) == '200':
            self.combobox_server.addItem('CYK')

        # buttons
        self.button_1 = QPushButton('Step 1: Retrieve Transcripts', self)
        self.button_1.resize(200, 20)
        self.button_1.move(20, 125)

        self.button_2 = QPushButton('Step 2: Start AutoPrimer', self)
        self.button_2.resize(200, 20)
        self.button_2.move(20, 150)
        self.button_2.setDisabled(True)

        # connect buttons to on_click functions
        self.button_1.clicked.connect(self.on_click_1)
        self.button_2.clicked.connect(self.on_click_2)

        # status bar
        self.statusBar().showMessage('Ready.')
        self.show()

        # output file
        self.outfile = ''

        # space-holder variables
        self.default_server = None
        self.gene = None
        self.transcript = None
        self.exon_regions = None
        self.translated_regions = None
        self.translated_region_names = None
        self.c_translated_regions = None
        self.c_translated_region_names = None
        self.pp = autoprimer.PrimerPool()

        # other settings
        self.cached_SNPfree = True
        ''' Auto-retry settings '''
        self.auto_retry = True
        self.auto_retry_level = 3
        ''' SNP and repeat masking settings '''
        self.routine_SNP_threshold = 0.0001
        self.retry_SNP_threshold = 'COMMON'
        self.routine_repeat_masking = 'hard'
        self.retry_repeat_masking = 'soft'
        ''' Amplicon size, exon-cutting and exon-combining settings '''
        self.routine_total_flank = 500
        self.retry_total_flank = 1000
        self.routine_target_flank = 100
        self.retry_target_flank = 80
        self.routine_amplicon_max = 600
        self.retry_amplicon_max = 800
        self.auto_cut_chunk = 400
        self.auto_combine_chunk = 500

    def on_click_1(self):
        self.statusBar().showMessage('Now getting transcript info...')
        gene_name = str(self.textbox.text()).upper()
        self.gene = autoprimer.Gene(gene_name)
        print(self.gene, file=sys.stderr)
        transcripts = self.gene.list_transcripts()
        print(transcripts, file=sys.stderr)
        QMessageBox.question(self, 'Transcript information', "The following transcripts are available: " + ', '.join(transcripts), QMessageBox.Ok, QMessageBox.Ok)
        self.textbox.setDisabled(True)
        self.button_1.setDisabled(True)
        for transcript in transcripts:
            self.combobox.addItem(transcript)
        self.combobox.setEnabled(True)
        self.button_2.setEnabled(True)
        self.statusBar().showMessage('ENSEMBL data retrieved.')


    def on_click_2(self):
        self.button_2.setDisabled(True)
        self.transcript = self.combobox.currentText()
        self.outfile = self.saveFileDialog(self.transcript + '.autoprimer.txt')
        self.default_server = self.combobox_server.currentText()
        ##################################################
        # START of AutoPrimer core logic
        ##################################################
        self.statusBar().showMessage('Now setting trasnscript as ' + self.transcript)
        QApplication.processEvents()
        self.gene.set_transcript(self.transcript)
        self.statusBar().showMessage('Now retrieving exon information...')
        QApplication.processEvents()
        self.exon_regions = self.gene.list_exon_regions()
        self.statusBar().showMessage('Now trimming exons to coding regions...')
        QApplication.processEvents()
        self.translated_regions = [self.gene.exon_to_translated(er) for er in self.exon_regions]
        self.translated_region_names = ['E' + str(i+1) for i in range(len(self.translated_regions))]
        self.statusBar().showMessage('Now splitting large amplicons...')
        QApplication.processEvents()
        self.translated_regions, self.translated_region_names = autoprimer.split_big_regions(self.translated_regions, self.translated_region_names, chunk=self.auto_cut_chunk)
        self.statusBar().showMessage('Now combining small amplicons...')
        QApplication.processEvents()
        reduction = 1
        iteration = 1
        while reduction > 0:
            self.statusBar().showMessage('Now combining small amplicons...(iteration ' + str(iteration) + ')' )
            QApplication.processEvents()
            iteration += 1
            self.c_translated_regions, self.c_translated_region_names = autoprimer.combine_adj_regions(self.translated_regions, self.translated_region_names, chunk=self.auto_combine_chunk)
            reduction = len(self.translated_regions) - len(self.c_translated_regions)
            if reduction > 0:
                self.translated_regions, self.translated_region_names = self.c_translated_regions, self.c_translated_region_names
            if reduction < 0:
                raise RuntimeError('An error occured. The number of combined regions is greated than input. Aborting.')
        self.statusBar().showMessage('Finished combining amplicons.')
        QApplication.processEvents()

        with open(self.outfile, 'w') as f:
            for translated_region, translated_region_name in zip(self.translated_regions, self.translated_region_names):
                if translated_region[1] and translated_region[2]:
                    self.statusBar().showMessage('Retrieving masked sequences... (' + translated_region_name + ')')
                    QApplication.processEvents()
                    upstream, downstream = autoprimer.get_flanking_regions(translated_region,
                                                                flank=self.routine_total_flank)
                    upseq = autoprimer.get_masked_sequence(upstream,
                                                self.default_server,
                                                threshold=self.routine_SNP_threshold,
                                                repeat_mask=self.routine_repeat_masking)
                    coreseq = autoprimer.get_sequence(translated_region)
                    downseq = autoprimer.get_masked_sequence(downstream,
                                                  self.default_server,
                                                  threshold=self.routine_SNP_threshold,
                                                  repeat_mask=self.routine_repeat_masking)
                    # prepare wrimer3 input
                    self.statusBar().showMessage('Generating primers... (' + translated_region_name + ')')
                    QApplication.processEvents()
                    sequence_id = self.gene.name + '_' + translated_region_name
                    sequence = upseq + coreseq + downseq
                    target_start = len(upseq) - self.routine_target_flank
                    target_length = len(coreseq) + 2 * self.routine_target_flank
                    target = str(target_start) + ',' + str(target_length)
                    wrimer3_result = autoprimer.wrimer3(sequence_id, sequence, target, self.default_server, max_prod=self.routine_amplicon_max)
                    primers = autoprimer.wrimer3_to_primers(wrimer3_result)
                    # if no primers found, retry using level 1 settings
                    if len(primers) == 0 and self.auto_retry == True and self.auto_retry_level >= 1:
                        self.statusBar().showMessage('Generating primers with relaxed settings level-1... (' + translated_region_name + ')')
                        QApplication.processEvents()
                        wrimer3_result = autoprimer.wrimer3(sequence_id, sequence, target, self.default_server, max_prod=self.retry_amplicon_max)
                        primers = autoprimer.wrimer3_to_primers(wrimer3_result)
                        # if no primers found, retry using level 2 settings
                        if len(primers) == 0 and self.auto_retry == True and self.auto_retry_level >= 2:
                            self.statusBar().showMessage('Generating primers with relaxed settings level-2... (' + translated_region_name + ')')
                            QApplication.processEvents()
                            upstream, downstream = autoprimer.get_flanking_regions(translated_region,
                                                            flank=self.retry_total_flank)
                            upseq = autoprimer.get_masked_sequence(upstream,
                                                        self.default_server,
                                                        threshold=self.routine_SNP_threshold,
                                                        repeat_mask=self.routine_repeat_masking)
                            coreseq = autoprimer.get_sequence(translated_region)
                            downseq = autoprimer.get_masked_sequence(downstream,
                                                          self.default_server,
                                                          threshold=self.routine_SNP_threshold,
                                                          repeat_mask=self.routine_repeat_masking)
                            sequence = upseq + coreseq + downseq
                            target_start = len(upseq) - self.retry_target_flank
                            target_length = len(coreseq) + 2 * self.retry_target_flank
                            target = str(target_start) + ',' + str(target_length)
                            wrimer3_result = autoprimer.wrimer3(sequence_id, sequence, target, self.default_server, max_prod=self.retry_amplicon_max)
                            primers = autoprimer.wrimer3_to_primers(wrimer3_result)
                            # if no primers found, retry using level 3 settings
                            if len(primers) == 0 and self.auto_retry == True and self.auto_retry_level >= 3:
                                self.statusBar().showMessage('Generating primers with relaxed settings level-3... (' + translated_region_name + ')')
                                QApplication.processEvents()
                                upstream, downstream = autoprimer.get_flanking_regions(translated_region,
                                                                flank=self.retry_total_flank)
                                upseq = autoprimer.get_masked_sequence(upstream,
                                                            self.default_server,
                                                            threshold=self.retry_SNP_threshold,
                                                            repeat_mask=self.retry_repeat_masking)
                                coreseq = autoprimer.get_sequence(translated_region)
                                downseq = autoprimer.get_masked_sequence(downstream,
                                                              self.default_server,
                                                              threshold=self.retry_SNP_threshold,
                                                              repeat_mask=self.retry_repeat_masking)
                                sequence = upseq + coreseq + downseq
                                target_start = len(upseq) - self.retry_target_flank
                                target_length = len(coreseq) + 2 * self.retry_target_flank
                                target = str(target_start) + ',' + str(target_length)
                                wrimer3_result = autoprimer.wrimer3(sequence_id, sequence, target, self.default_server, max_prod=self.retry_amplicon_max)
                                primers = autoprimer.wrimer3_to_primers(wrimer3_result)
                                if len(primers) == 0:
                                    print('# ' + translated_region_name, file=f)
                                    print('# No primers found after 3 levels of intelligent-retry. Proceeding to next exon...', file=f)
                    if self.cached_SNPfree:
                        primer_count = 1
                        pair_name, fp, rp = list(), list(), list()
                        for p in primers:
                            pair_name.append(sequence_id + '_' + str(primer_count))
                            fp.append(p[0])
                            rp.append(p[1])
                        self.statusBar().showMessage('SNPfree analysis... (initial submission)')
                        QApplication.processEvents()
                        autoprimer.autoSNPfree(pair_name, fp, rp, self.default_server, batch_mode=True)

                    primer_count = 1
                    for p in primers:
                        pair_name = sequence_id + '_' + str(primer_count)
                        fp = p[0]
                        rp = p[1]
                        print(pair_name, fp, rp, sep='\t', file=f)
                        self.statusBar().showMessage('SNPfree analysis... (' + pair_name + ')')
                        QApplication.processEvents()
                        autoSNPfree_result = autoprimer.autoSNPfree(pair_name, fp, rp, self.default_server)
                        primer_name, primer_score = autoprimer.SNPfree_to_score(autoSNPfree_result)
                        print('#', primer_name, 'SNPfree score:', primer_score, file=f)
                        self.statusBar().showMessage('SNPfree analysis... (' + pair_name + ')' + ' Score: ' + str(primer_score))
                        QApplication.processEvents()
                        self.pp.add(sequence_id, primer_name, fp, rp,
                                    primer_score, autoprimer.SNPfree_to_score(autoSNPfree_result, get_size=True))
                        primer_count += 1
                else:
                    self.statusBar().showMessage('Skipping non-translated exon... (' + translated_region_name + ')')
                    QApplication.processEvents()

            self.statusBar().showMessage('Exporting results...')
            QApplication.processEvents()
            print('##############################')
            print('# Below is the list of selected primers for the target regions:', file=f)
            final_list = self.pp.getbest()
            for pair in final_list:
                print('#', pair[1], 'chosen for', pair[0], 'with SNPfree score', pair[4], file=f)
                print(pair[0], pair[2], pair[3], sep='\t', file=f)

            print('##############################', file=f)
            print('# Below is the list for ordering:', file=f)
            for pair in final_list:
                p_gene, p_region = str(pair[0]).split('_')
                print(pair[2], len(pair[2]), 'No', p_gene, p_region, 'F', pair[5], sep='\t', file=f)
                print(pair[3], len(pair[3]), 'No', p_gene, p_region, 'R', pair[5], sep='\t', file=f)
            self.statusBar().showMessage('Finished!')
            QApplication.processEvents()
        ##################################################
        # END of AutoPrimer core logic
        ##################################################
    def saveFileDialog(self, default_filename):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Choose AutoPrimer output location", default_filename,"All Files (*);;Text Files (*.txt)", options=options)
        if fileName:
            return fileName
        else:
            self.on_click_1()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    app.exec_()