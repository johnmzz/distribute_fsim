# Importing difflib
import difflib
values = {}

with open('result.txt') as file_1:
    lines = file_1.readlines()
    for line in lines:
        words = line.strip().split(" ")
        if (words[0],words[1]) not in values.keys():
            values[(words[0],words[1])] = words[2]


with open('result_partition1.txt') as file_2:
    lines = file_2.readlines()
    for line in lines:
        words = line.strip().split(" ")
        if (words[0], words[1]) in values.keys():
            assert words[2] == values[(words[0],words[1])