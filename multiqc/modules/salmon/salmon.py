#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os
from .GCModel import GCModel
from .SeqModel import SeqModel
from multiqc.plots import linegraph
from multiqc.plots import heatmap
from multiqc.modules.base_module import BaseMultiqcModule
import numpy as np

# Function to check if the directory contains GC or Seq Bias files
def checkJSONForBias(directory, checkBias):
    is_file = False
    for fname in os.listdir(directory):
        if fname.endswith('.json'):
            filename = os.path.sep.join([directory, fname])
            with open(filename, 'r') as f:
                jsonContents = json.load(f)
                if checkBias in jsonContents:
                    is_file = True
                    return is_file

    return is_file


# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Salmon', anchor='salmon',
        href='http://combine-lab.github.io/salmon/',
        info="is a tool for quantifying the expression of transcripts using RNA-seq data.")

        # Parse meta information. JSON win!
        self.salmon_meta = dict()

        # Declaring dicts to hold ratios for first,midddle,last rows with weights and avergage ratio for GC Bias
        self.salmon_bias_FirstSampleWeights = dict()
        self.salmon_bias_MiddleSampleWeights = dict()
        self.salmon_bias_LastSampleWights = dict()
        self.salmon_bias_Average=dict()
        self.salmon_bias_TotalAverage = dict()
        self.salmon_bias_Merged = []
        self.heatmap_gc_Merged = []
        self.heatmap_seq_Merged = []

        #Declaring dicts to hold sequence 3' and 5' marginalized ratio for all bases i.e A,C,G,T and the average bias for 3' and 5'
        self.salmon_seq3A = dict()
        self.salmon_seq3C = dict()
        self.salmon_seq3G = dict()
        self.salmon_seq3T = dict()
        self.salmon_seq3Average = dict()
        self.salmon_seq3_Merged = []
        self.salmon_seq5A = dict()
        self.salmon_seq5C = dict()
        self.salmon_seq5G = dict()
        self.salmon_seq5T = dict()
        self.salmon_seq5Average = dict()
        self.salmon_seq5_Merged = []

        # Declaring lists to hold arrays for every sample used in Heatmaps
        self.heatmapFirstrow=[]
        self.heatMapMiddleRow=[]
        self.heatMapLastRow=[]
        self.averageBiasHeatMap=[]
        self.salmon_seq3HeatMap=[]
        self.salmon_seq5HeatMap=[]

        # List of all the sample names
        self.sample_names=[]

        count = 0
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename( os.path.dirname(f['root']) )
            s_name = self.clean_s_name(s_name, f['root'])
            self.salmon_meta[s_name] = json.loads(f['f'])
            s_name_trimmed = s_name.partition('|')[0].split()
            self.sample_names.append(s_name_trimmed)

            # Check if folder contains GC bias files
            gcBias = checkJSONForBias(os.path.dirname(f['root']), 'gcBias')

            if gcBias:
                # Dicts for every sample for all the bucket(25) ratios to hold (x,y) data for linegraphs
                firstRatioWeight = OrderedDict()
                middleRatioWeight = OrderedDict()
                lastRatioWeight = OrderedDict()
                average = OrderedDict()
                sampleAverage = OrderedDict()

                gc = GCModel() # Instantiate GCModel class
                # Call the GCModel method to get all observed and expected values
                gc.from_file(os.path.dirname(f['root']))
                first_Row = (gc.obs_[0] / gc.exp_[0])*(gc.obs_weights_[0]/gc.exp_weights_[0])
                middle_Row = (gc.obs_[1] / gc.exp_[1])*(gc.obs_weights_[1]/gc.exp_weights_[1])
                last_Row = (gc.obs_[2] / gc.exp_[2])*(gc.obs_weights_[2]/gc.exp_weights_[2])
                #Avergaing ratios for each row used in Heatmap for every row
                self.heatmapFirstrow.append(first_Row.tolist())
                self.heatMapMiddleRow.append(middle_Row.tolist())
                self.heatMapLastRow.append(last_Row.tolist())

                heatmapAverage=[]
                # Iterating over all the buckets to create Ordered Dicts
                for i in range(len(first_Row)):
                    index = i*(100/len(first_Row))
                    firstRatioWeight[index] = first_Row[i]
                    middleRatioWeight[index] = middle_Row[i]
                    lastRatioWeight[index] = last_Row[i]
                    average[index]=np.mean([first_Row[i],middle_Row[i],last_Row[i]])
                    heatmapAverage.append(average[index])

                # Setting all the ordered dicts to the outermost Dictionaries with sample name as keys
                self.salmon_bias_FirstSampleWeights[s_name]=firstRatioWeight
                self.salmon_bias_MiddleSampleWeights[s_name] = middleRatioWeight
                self.salmon_bias_LastSampleWights[s_name] = lastRatioWeight
                self.salmon_bias_Average[s_name] = average
                self.averageBiasHeatMap.append(heatmapAverage)

            # Check if folder contains sequence bias files
            seqBias = checkJSONForBias(os.path.dirname(f['root']), 'seqBias')
            if seqBias:
                # Dicts for every base for 3' and 5' sequence, average 3' and average 5' and quant dict
                seq3A = OrderedDict()
                seq5A = OrderedDict()
                seq3C = OrderedDict()
                seq5C = OrderedDict()
                seq3G = OrderedDict()
                seq5G = OrderedDict()
                seq3T = OrderedDict()
                seq5T = OrderedDict()
                seq3_Average = OrderedDict()
                seq5_Average = OrderedDict()

                # Calculate the ratio of all rows for observed by expected
                seq = SeqModel()# Instantiate SeqModel class
                # Call the SeqModel method to get all observed and expected ratios
                seq.from_file(os.path.dirname(f['root']))
                seq3A_prob = seq.obs3_[0] / seq.exp3_[0]
                seq3C_prob = seq.obs3_[1] / seq.exp3_[1]
                seq3G_prob = seq.obs3_[2] / seq.exp3_[2]
                seq3T_prob = seq.obs3_[3] / seq.exp3_[3]
                seq5A_prob = seq.obs5_[0] / seq.exp5_[0]
                seq5C_prob = seq.obs5_[1] / seq.exp5_[1]
                seq5G_prob = seq.obs5_[2] / seq.exp5_[2]
                seq5T_prob = seq.obs5_[3] / seq.exp5_[3]

                seq3_HeatMap = []
                seq5_HeatMap = []
                # Iterate over the contect length to create all Orderede Dictonaries of (x,y) values for linegraph and list for Heatmap
                for i in range(len(seq3A_prob)):
                    index = i * (100 / len(seq3A_prob))
                    seq3A[index] = seq3A_prob[i]
                    seq5A[index] = seq5A_prob[i]
                    seq3C[index] = seq3C_prob[i]
                    seq5C[index] = seq5C_prob[i]
                    seq3G[index] = seq3G_prob[i]
                    seq5G[index] = seq5G_prob[i]
                    seq3T[index] = seq3T_prob[i]
                    seq5T[index] = seq5T_prob[i]
                    seq3_Average[index] = np.mean([seq3A_prob[i], seq3C_prob[i], seq3G_prob[i], seq3T_prob[i]])
                    seq5_Average[index] = np.mean([seq5A_prob[i], seq5C_prob[i], seq5G_prob[i], seq5T_prob[i]])
                    seq3_HeatMap.append(seq3_Average[index])
                    seq5_HeatMap.append(seq5_Average[index])

                # Setting all the ordered dicts to the outermost Dictionaries with sample name as keys
                self.salmon_seq3A[s_name] = seq3A
                self.salmon_seq5A[s_name] = seq5A
                self.salmon_seq3C[s_name] = seq3C
                self.salmon_seq5C[s_name] = seq5C
                self.salmon_seq3G[s_name] = seq3G
                self.salmon_seq5G[s_name] = seq5G
                self.salmon_seq3T[s_name] = seq3T
                self.salmon_seq5T[s_name] = seq5T
                self.salmon_seq3Average[s_name] = seq3_Average
                self.salmon_seq5Average[s_name] = seq5_Average
                self.salmon_seq3HeatMap.append(seq3_HeatMap)
                self.salmon_seq5HeatMap.append(seq5_HeatMap)

        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        for f in self.find_log_files('salmon/fld'):
            # Get the s_name from the parent directory
            if os.path.basename(f['root']) == 'libParams':
                s_name = os.path.basename( os.path.dirname(f['root']) )
                s_name = self.clean_s_name(s_name, f['root'])
                parsed = OrderedDict()
                for i, v in enumerate(f['f'].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed

        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning

        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            self.write_data_file(self.salmon_meta, 'multiqc_salmon')
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        if len(self.salmon_bias_Average) > 0:
            log.info("Found {} GC Bias".format(len(self.salmon_bias_Average)))

        if len(self.salmon_seq3Average) > 0:
            log.info("Found {} Sequence 3' bias".format(len(self.salmon_seq3Average)))

        if len(self.salmon_seq5Average) > 0:
            log.info("Found {} Sequence 5' bias".format(len(self.salmon_seq5Average)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['percent_mapped'] = {
            'title': '% Aligned',
            'description': '% Mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['num_mapped'] = {
            'title': 'M Aligned',
            'description': 'Mapped reads (millions)',
            'min': 0,
            'scale': 'PuRd',
            'modify': lambda x: float(x) / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.salmon_meta, headers)

        # Fragment length distribution plot
        pconfig = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Fragment Length Distribution',
            'ylab': 'Fraction',
            'xlab': 'Fragment Length (bp)',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        self.add_section( plot = linegraph.plot(self.salmon_fld, pconfig) )

        # GC Bias Plots.

        filenames = self.salmon_bias_FirstSampleWeights.keys()
        rows = []
        count = 0
        plot = {}
        data_labels=[]
        for key in filenames:
            temp = {}
            for keys in self.salmon_bias_FirstSampleWeights[key]:
                temp[keys] = self.salmon_bias_FirstSampleWeights[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_bias_Merged.append(plot)
        data_labels.append({'name':'first'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_bias_MiddleSampleWeights[key]:
                temp[keys] = self.salmon_bias_MiddleSampleWeights[key][keys]
            rows.append(temp)    
            plot[key] = rows[count]
            count += 1
        self.salmon_bias_Merged.append(plot)
        data_labels.append({'name':'middle'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_bias_LastSampleWights[key]:
                temp[keys] = self.salmon_bias_LastSampleWights[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_bias_Merged.append(plot)
        data_labels.append({'name':'last'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_bias_Average[key]:
                temp[keys] = self.salmon_bias_Average[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_bias_Merged.append(plot)
        data_labels.append({'name':'average'})

        pconfig_bias_Merged = {
            'smooth_points': 500,
            'title': 'GC Bias Plots.',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'xmax': 100,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
            'data_labels': data_labels
        }
        self.add_section(name='GC Bias Plots',plot=linegraph.plot(self.salmon_bias_Merged, pconfig_bias_Merged))

        # Seq 3 Plots.

        filenames = self.salmon_seq3A.keys()
        rows = []
        count = 0
        plot = {}
        data_labels=[]
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq3A[key]:
                temp[keys] = self.salmon_seq3A[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq3_Merged.append(plot)
        data_labels.append({'name':'A'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq3C[key]:
                temp[keys] = self.salmon_seq3C[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq3_Merged.append(plot)
        data_labels.append({'name':'C'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq3G[key]:
                temp[keys] = self.salmon_seq3G[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq3_Merged.append(plot)
        data_labels.append({'name':'G'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq3T[key]:
                temp[keys] = self.salmon_seq3T[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq3_Merged.append(plot)
        data_labels.append({'name':'T'})

        rows = []
        temp={}
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq3Average[key]:
                temp[keys] = self.salmon_seq3Average[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq3_Merged.append(plot)
        data_labels.append({'name':'Average'})

        pconfig_seq3_Merged = {
            'smooth_points': 500,
            'title': 'Seq 3 Plots.',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'xmax': 100,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
            'data_labels': data_labels
        }
        self.add_section(name='Seq 3 Plots',plot=linegraph.plot(self.salmon_seq3_Merged, pconfig_seq3_Merged))
        
        # Seq 5 Plots.

        filenames = self.salmon_seq5A.keys()
        rows = []
        count = 0
        plot = {}
        data_labels=[]
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq5A[key]:
                temp[keys] = self.salmon_seq5A[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq5_Merged.append(plot)
        data_labels.append({'name':'A'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq5C[key]:
                temp[keys] = self.salmon_seq5C[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq5_Merged.append(plot)
        data_labels.append({'name':'C'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq5G[key]:
                temp[keys] = self.salmon_seq5G[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq5_Merged.append(plot)
        data_labels.append({'name':'G'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq5T[key]:
                temp[keys] = self.salmon_seq5T[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq5_Merged.append(plot)
        data_labels.append({'name':'T'})

        rows = []
        count = 0
        plot = {}
        for key in filenames:
            temp = {}
            for keys in self.salmon_seq5Average[key]:
                temp[keys] = self.salmon_seq5Average[key][keys]
            rows.append(temp)
            plot[key] = rows[count]
            count += 1
        self.salmon_seq5_Merged.append(plot)
        data_labels.append({'name':'Average'})

        pconfig_seq5_Merged = {
            'smooth_points': 500,
            'title': 'Seq 5 Plots.',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'xmax': 100,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
            'data_labels': data_labels
        }
        self.add_section(name='Seq 5 Plots',plot=linegraph.plot(self.salmon_seq5_Merged, pconfig_seq5_Merged))

        # First Row of all samples Heatmap
        FirstRowCoff = np.corrcoef(self.heatmapFirstrow)
        self.add_section(name='GC Bias First Row Heatmap', description='Heatmap to display variance between first row ratios of all the samples',
                         plot=heatmap.plot(FirstRowCoff, self.sample_names, self.sample_names))

        # Middle Row of all samples Heatmap
        MiddleRowCoff = np.corrcoef(self.heatMapMiddleRow)
        self.add_section(name='GC Bias Middle Row Heatmap', description='Heatmap to display variance between middle row ratios of all the samples',
                         plot=heatmap.plot(MiddleRowCoff, self.sample_names, self.sample_names))

        # Last Row of all samples Heatmap
        LastRowCoff = np.corrcoef(self.heatMapLastRow)
        self.add_section(name='GC Bias Last Row Heatmap', description='Heatmap to display variance between last row ratios of all the samples',
                         plot=heatmap.plot(LastRowCoff, self.sample_names, self.sample_names))

        # GC Bias HeatMap
        AverageCoff = np.corrcoef(self.averageBiasHeatMap)
        self.add_section(name='GC Bias Heatmap', description='Heatmap to display average bias across all samples',
                         plot=heatmap.plot(AverageCoff, self.sample_names, self.sample_names))

        # Seq 3' Heatmap
        Seq3HeatMap = np.corrcoef(self.salmon_seq3HeatMap)
        self.add_section(name='Sequence 3 Heatmap', description='Heatmap to display Sequence 3 prime across all samples',
                         plot=heatmap.plot(Seq3HeatMap, self.sample_names, self.sample_names))

        # Seq 5' Heatmap
        Seq5HeatMap = np.corrcoef(self.salmon_seq5HeatMap)
        self.add_section(name='Sequence 5 Heatmap', description='Heatmap to display Sequence 5 prime across all samples',
                         plot=heatmap.plot(Seq5HeatMap, self.sample_names, self.sample_names))
