import numpy as np

class SeqModel:
    def __init__(self):
        self.obs3_ = None
        self.exp3_ = None
        self.obs5_ = None
        self.exp5_ = None
        self.valid_ = False

    def populate_model_(self, data_):
        import struct
        import time

        offset = 0
        int_struct = struct.Struct('@i')
        long_struct = struct.Struct('@q')

        context_length = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        second = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        third = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        firstArray_struct = struct.Struct('@' + context_length*'i')
        first_Array = firstArray_struct.unpack_from(data_[offset:])
        offset += firstArray_struct.size

        secondArray_struct = struct.Struct('@' + context_length*'i')
        second_Array = secondArray_struct.unpack_from(data_[offset:])
        offset += secondArray_struct.size

        thirdArray_struct = struct.Struct('@' + context_length*'i')
        third_Array = thirdArray_struct.unpack_from(data_[offset:])
        offset += thirdArray_struct.size

        nrow = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        ncol = long_struct.unpack_from(data_[offset:])[0]
        offset += long_struct.size

        vlmm_struct = struct.Struct('@' + nrow * ncol * 'd')
        vlmm = vlmm_struct.unpack_from(data_[offset:])
        vlmm = np.array(vlmm)
        vlmm = vlmm.reshape(ncol, nrow).T
        vlmm = (vlmm.T / vlmm.sum(axis=1)).T
        offset += vlmm_struct.size

        nrow_1 = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        ncol_1 = int_struct.unpack_from(data_[offset:])[0]
        offset += int_struct.size

        margin_struct = struct.Struct('@' + nrow_1 * context_length * 'd')
        margin = margin_struct.unpack_from(data_[offset:])

        margin = np.array(margin)
        margin = margin.reshape(context_length, nrow_1).T
        margin = (margin.T / margin.sum(axis=1)).T

        return margin


    # dname is the root directory of salmon output
    def from_file(self, dname):
        import os
        import gzip

        obs3_name = os.path.sep.join([dname, 'aux_info', 'obs3_seq.gz'])
        exp3_name = os.path.sep.join([dname, 'aux_info', 'exp3_seq.gz'])

        obs5_name = os.path.sep.join([dname, 'aux_info', 'obs5_seq.gz'])
        exp5_name = os.path.sep.join([dname, 'aux_info', 'exp5_seq.gz'])

        # Observed and Expected for Sequence 3' Bias ############################
        obs3_dat, exp3_dat = None, None
        try:
            with gzip.open(obs3_name) as obs3_file:
                obs3_dat = obs3_file.read()
            self.obs3_ = self.populate_model_(obs3_dat)
            #print(self.obs3_)
        except IOError:
            print("Could not open file {}".format(obs3_name))
            return False

        try:
            with gzip.open(exp3_name) as exp3_file:
                exp3_dat = exp3_file.read()
            self.exp3_ = self.populate_model_(exp3_dat)
        except IOError:
            print("Could not open file {}".format(exp3_name))
            return False

        # Observed and Expected for Sequence 5' Bias ############################
        obs5_dat, exp5_dat = None, None
        try:
            with gzip.open(obs5_name) as obs5_file:
                obs5_dat = obs5_file.read()
            self.obs5_ = self.populate_model_(obs5_dat)
        except IOError:
            print("Could not open file {}".format(obs5_name))
            return False

        try:
            with gzip.open(exp5_name) as exp5_file:
                exp5_dat = exp5_file.read()
            self.exp5_ = self.populate_model_(exp5_dat)
        except IOError:
            print("Could not open file {}".format(exp5_name))
            return False


        self.valid_ = True
        return True
