#!/usr/bin/env python
# Modified propagate.py to work with SALSA Workflow Manager

import sys
import argparse
import numpy as np
import pandas as pd
import re
import os
import shutil
import subprocess
import base64
import datetime
import collections

timestamp_format = "%H:%M:%S %Y-%m-%d"

def generate_unique_id(length=5):
    return base64.b64encode(os.urandom(64)).decode().replace("/", "").replace("+", "")[:length]

def get_time():
    return datetime.datetime.now().strftime(timestamp_format)

def get_file_timestamp(path):
    int_timestamp = os.path.getmtime(path)
    return datetime.datetime.fromtimestamp(int_timestamp).strftime(timestamp_format)

def str_to_dict(dict_string, value_type="int"):
    new_list = dict_string.split("::")
    new_dict = {}
    for pair in new_list:
        if ":" in pair:
            key, value = pair.split(":")
            if value_type == "int":
                try:
                    value = int(value)
                except ValueError:
                    pass
            new_dict[key] = value
    return new_dict

def dict_to_str(a_dict):
    dict_string = ""
    for key, value in a_dict.items():
        dict_string += "{}:{}::".format(key, value)
    return dict_string[:-2]

class Logger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.error = sys.stderr
        self.logfile = open(filename, "a")
        sys.stdout = self
        sys.stderr = self

    def write(self, message, logfile_only=False):
        if not logfile_only:
            self.terminal.write(message)
        self.logfile.write(message)

    def flush(self):
        pass

    def stop(self):
        sys.stdout.logfile.close()
        sys.stdout = sys.stdout.terminal
        sys.stderr = sys.stderr.error

class InventoryRow:
    def __init__(self, pd_series):
        self.import_attributes_from_existing_row(pd_series)
