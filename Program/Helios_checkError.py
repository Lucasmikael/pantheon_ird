import os
import csv

def checkFormatGeneList(filename):
    check = True
    check_len = True
    check_row = True
    check_error = True
    file = open(filename, "r")
    reader = csv.reader(file, delimiter="\t")
    ncol = len(next(reader))
    for row in reader :
        if ncol < 3 :
            check_len = False
        else:
            first_row = isinstance(row[0], str)
            second_row = isinstance(row[1], str)
            try :
                third_row = int(row[2])
                if first_row == False or second_row == False:
                    check_row = False
            except :
                check = False
        if check_len == False or check_row == False or check_error == False :
            check = False
            break
    return check


def checkFormatInteraction(filename):
    check = True
    check_len = True
    check_int = False
    check_row = True
    check_error = True
    file = open(filename, "r")
    reader = csv.reader(file, delimiter="\t")
    ncol = len(next(reader))
    next(reader)
    for row in reader :
        if ncol < 6 :
            check_len = False
        else :
            first_row = isinstance(row[0], str)
            try :
                second_row = int(row[1])
                third_row = isinstance(row[2], str)
                if first_row == False or third_row == False :
                    check_row = False
                if second_row == 1 or second_row == -1:
                    check_int = True
            except :
                check_error = False
        if check_len == False or check_int == False or check_row == False or check_error == False :
            check = False
            break
    return check

def checkLoadElement(filename) :
    check = True
    check_len = True
    check_int = False
    check_row = True
    check_error = True
    file = open(filename, "r")
    reader = csv.reader(file, delimiter=",")
    next(reader)
    ncol = len(next(reader))
    for row in reader :
        if ncol < 3 :
            check_len = False
        else:
            first_row = isinstance(row[0], str)
            third_row = isinstance(row[2], str)
            try :
                second_row = int(row[1])
                if first_row == False or third_row == False:
                    check_row = False
                if second_row == 1 or second_row == -1:
                    check_int = True
            except :
                check_error = False
        if check_len == False or check_int == False or check_row == False or check_error == False :
            check = False
            break
    return check



def checkSavedNetwork(filename) :
    check = True
    check_len = True
    check_int = False
    check_row = True
    check_error = True
    file = open(filename, "r")
    reader = csv.reader(file, delimiter=",")
    next(reader)
    ncol = len(next(reader))
    for row in reader:
        if ncol < 4 :
            check_len = False
        if row[0] == "Etat" :
            break
        first_row = isinstance(row[0], str)
        second_row = isinstance(row[1], str)
        try :
            third_row = int(row[2])
            fourth_row = int(row[3])
            if first_row == False or second_row == False:
                check_row = False
            if third_row == 1 or third_row == -1 and fourth_row == 0 or fourth_row == 1 :
                check_int = True
        except :
            check_error = False
        if check_len == False or check_int == False or check_row == False or check_error == False :
            check = False
            break

    return check

def checkPath(filename):
    value_error_path = 0
    check = os.path.isfile(filename)
    if check == False :
        value_error_path = 2

    return check, value_error_path
