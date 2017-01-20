# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 21:34:24 2017

@author: Dell
"""

#!/bin/env python3
import sys
import struct
import numpy
import pprint

def parseCEL(filename):
    if filename == None:
        print("\nUSAGE\n\n\
    Pass the name of the CEL file as the first and only argument to the script. Execute the script with python3, for example:\n\n\
    python3 parseCELfile.py 329975487_M.CEL\n")
        return
    numbers = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
    numbersMap = {}
    headerMap = {}
    
    with open(filename, "rb") as f:
    
        for name in numbers:
            numbersMap[name] = struct.unpack('<i', f.read(4))[0]
        char = f.read(numbersMap["headerLen"])
        header = char.decode("ascii", "ignore")
        for header in header.split("\n"):
            if "=" in header:
                header = header.split("=")
                headerMap[header[0]] = header[1]
        
        char = b'\x00'
        safetyValve = 10**4
        for i in range(10**4):
            char = f.read(1)
            if char == b'\x04':
                break
            if i == safetyValve:
                raise(ValueError("Parse Error"))
                
        padding = f.read(15)
    
        structa = struct.Struct("< f f h")
    
        structSize = 10
    
        length = numbersMap["cellNo"]
    
        intensities = numpy.empty(length, dtype=float)
        deviations = numpy.empty(length, dtype=float)
        pixels = numpy.empty(length, dtype=int)
    
        b = f.read(structSize * length)
        for i in range(length):
            binaryFragment = b[i * structSize : (i + 1) * structSize]
            intensity, deviation, pixelcount = structa.unpack(binaryFragment)
            intensities[i] = intensity
            deviations[i] = deviation
            pixels[i] = pixelcount
    
    return (intensities, deviations, pixels)