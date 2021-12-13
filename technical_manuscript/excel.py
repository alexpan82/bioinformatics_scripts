import csv
import openpyxl
from sys import argv
import pandas as pd
from datetime import date
script= argv[0]
overlap_name = argv[1]
Patient = argv[2]

input_file = overlap_name
output_file = "Table for " + Patient +"_spliced_withoutUTRs.xlsx"

wb = openpyxl.load_workbook(output_file)
ws = wb.create_sheet("Overlap_list")


with open(input_file, 'rb') as data:
    
    for row in data:
        ws.append(row.split(','))

wb.save(output_file)