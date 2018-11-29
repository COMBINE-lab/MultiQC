#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
import json
import logging
import os
from .GCModel import GCModel

from multiqc.plots import linegraph
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

        #Declaring dicts to hold sequence 3' and 5' marginalized ratio for all bases i.e A,C,G,T and the average bias for 3' and 5'
        self.salmon_seq3A = dict()
        self.salmon_seq3C = dict()
        self.salmon_seq3G = dict()
        self.salmon_seq3T = dict()
        self.salmon_seq5A = dict()
        self.salmon_seq5C = dict()
        self.salmon_seq5G = dict()
        self.salmon_seq5T = dict()
        self.salmon_seq3Average = dict()
        self.salmon_seq5Average = dict()

        #Declaring dict to hold the ratios of Effective v/s Actual length of samples from quant.sf file
        self.salmon_quant =dict()

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
                #print(gc.obs_)
                first_Row = (gc.obs_[0] / gc.exp_[0])*(gc.obs_weights_[0]/gc.exp_weights_[0])
                middle_Row = (gc.obs_[1] / gc.exp_[1])*(gc.obs_weights_[1]/gc.exp_weights_[1])
                last_Row = (gc.obs_[2] / gc.exp_[2])*(gc.obs_weights_[2]/gc.exp_weights_[2])

                # Avergaing all the ratios for the entire sample
                totalSampleAverage = ( (sum(first_Row)+sum(middle_Row)+sum(last_Row))/(len(first_Row)+len(middle_Row)+len(last_Row)))
                sampleAverage[count]=totalSampleAverage
                count = count +1
                self.salmon_bias_TotalAverage[s_name_trimmed[0]] = sampleAverage
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
        keyarray = self.salmon_bias_FirstSampleWeights.keys()
        temp2 = []
        temp={}
        count = 0
        final = {}
        data_labels=[]
        for key in keyarray:
            temp = {}
            for keys in self.salmon_bias_FirstSampleWeights[key]:
                temp[keys] = self.salmon_bias_FirstSampleWeights[key][keys]
            temp2.append(temp)
            final[key] = temp2[count]
            count += 1
        self.salmon_bias_Merged.append(final)
        data_labels.append({'name':'first'})

        temp2 = []
        temp={}
        count = 0
        final = {}
        for key in keyarray:
            temp = {}
            for keys in self.salmon_bias_MiddleSampleWeights[key]:
                temp[keys] = self.salmon_bias_MiddleSampleWeights[key][keys]
            temp2.append(temp)    
            final[key] = temp2[count]
            count += 1
        self.salmon_bias_Merged.append(final)
        data_labels.append({'name':'middle'})

        temp2 = []
        temp={}
        count = 0
        final = {}
        for key in keyarray:
            temp = {}
            for keys in self.salmon_bias_LastSampleWights[key]:
                temp[keys] = self.salmon_bias_LastSampleWights[key][keys]
            temp2.append(temp)
            final[key] = temp2[count]
            count += 1
        self.salmon_bias_Merged.append(final)
        data_labels.append({'name':'last'})

        temp2 = []
        temp={}
        count = 0
        final = {}
        for key in keyarray:
            temp = {}
            for keys in self.salmon_bias_Average[key]:
                temp[keys] = self.salmon_bias_Average[key][keys]
            temp2.append(temp)
            final[key] = temp2[count]
            count += 1
        self.salmon_bias_Merged.append(final)
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
