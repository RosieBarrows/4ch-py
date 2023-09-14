#!/usr/bin/env python2

import os
import sys
import numpy as np
import copy

import common_4ch.file_utils as fu
from common_4ch.config import configure_logging
from common_4ch.process_handler import correct_fibres
milog = configure_logging(log_name=__name__)

def main(mshName) :
	correct_fibres(mshName)

if __name__ == '__main__':
	main(sys.argv[1])