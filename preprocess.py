# -*- coding: utf-8 -*-
# @Time    : 2021/4/15 下午6:26
# @Author  : xiaorui su
# @Email   :  suxiaorui19@mails.ucas.edu.cn
# @File    : preprocess.py
# @Software : PyCharm

##dataset: HDVD
import untangle
import pandas as pd
import numpy as np
import os
from yoctol_utils.hash import consistent_hash


HDVD_path = "raw_data/HDVD/dataset/"
drug_smile = "drugs.csv"
virusdrug = "virusdrug.csv"


def drug_SMILE_read(file_path):

    drug_smile = []  #drugs,cas number, a number, smile
    f = open(file_path, 'r')
    flag = 1
    for line in f.readlines():
        if flag == 1:
            flag = 0
            continue
        print(line)
        row = line.strip().split(",")
        print(row)
        drug_smile.append(row)

    drug_smile = np.array(drug_smile)
    print(drug_smile.shape)
    drug_hash = convert_hash(drug_smile)
    #print(drug_smile)

    return drug_hash;

def convert_hash(drug_smile):

    smiles = []
    hashes = []
    for i in drug_smile:
        strings = []
        for c in i[3]:
            strings.append(str(c))
        smiles.append(strings)
        hashes.append(hash_seq(strings, 512))

    print(len(smiles))
    print(len(hashes))

    new_hash = []
    for i in hashes:
        if len(i) == 512:
            print("yes")
            new_hash.append(i)
        else:
            if len(i) < 512:
                padding = np.zeros(512 - len(i))
                new_hashes = np.concatenate((padding, i), axis=0)
                new_hash.append(new_hashes)
            else:
                new_hash.append(i[0:512])

    new_hash = np.array(new_hash)
    print(np.array(new_hash).shape)
    print(new_hash.reshape(len(hashes), 512))

    np.save("data/HDVD_smile_hash", new_hash)
    np.save("data/HDVD_id_smile", smiles)

    return new_hash;

def hash_seq(sequence, max_index):
    print(sequence)
    return np.array([consistent_hash(word) % max_index + 1 for word in sequence])

def load_interactions(aj_path):

    interactions = [] #drug id, virus id
    drug_name = []
    negative_interactions = []

    entity = "entity2id.txt"
    entity_path = "raw_data/HDVD/entity2id.txt"
    train_path = "raw_data/HDVD/train2id.txt"
    approved_path = "raw_data/HDVD/approved_example.txt"

    f = open(aj_path, 'r')
    flag = 1
    count = 0
    for line in f.readlines():
        if flag == 1:
            virus_name = line.strip().split(",")[1:]
            print(virus_name)
            print(len(virus_name))
            flag = 0
            continue

        row = line.strip().split(",")
        drug_name.append(row[0])
        count = count + 1
        for i in range(0,len(row)):
            if row[i] == '1':
                interactions.append([count,i+219])
            else:
                negative_interactions.append([count,i+219])

    interactions = np.array(interactions)
    negative_interactions = np.array(negative_interactions)
    print(interactions.shape)
    write_entity(entity_path,drug_name,virus_name)
    #write_train2id(train_path,interactions)
    #write_approved_example(approved_path, interactions,negative_interactions)

    return;

def write_entity(entity_path,drug_name,virus_name):

    for i in virus_name:
        print(i)
        drug_name.append(i)
    entity_name = drug_name
    print(entity_name)
    f = open(entity_path, 'a')
    f.write(str(len(entity_name)) + "\n")
    for i in range(0, len(entity_name)):
        f.write(str(entity_name[i]) + " ")
        f.write(str(i+1))
        f.write("\n")
    f.close()

    return;

def write_train2id(train_path,interactions):

    f = open(train_path,'a')
    f.write(str(interactions.shape[0]) + "\n")
    for i in range(0, interactions.shape[0]):
        f.write(str(interactions[i][0]) + " ")
        f.write(str(interactions[i][1]) + " ")
        f.write(str(1))
        f.write("\n")
    f.close()

    return

def write_approved_example(approved_path,interactions, negative_interactions):


    ##generate negative sampling
    drug_number = [] #1-219
    virus_number = [] #220-253

    rand_arr = np.arange(negative_interactions.shape[0])

    np.random.shuffle(rand_arr)

    generating = negative_interactions[rand_arr[0:interactions.shape[0]]]

    print(negative_interactions[rand_arr[0:interactions.shape[0]]])

    f = open(approved_path,'a')
    for i in range(0, interactions.shape[0]):
        f.write(str(interactions[i][0]) + " ")
        f.write(str(interactions[i][1]) + " ")
        f.write(str(1))
        f.write("\n")
    for i in range(0, generating.shape[0]):
        f.write(str(generating[i][0]) + " ")
        f.write(str(generating[i][1]) + " ")
        f.write(str(0))
        f.write("\n")
    f.close()

    return

def process_virus_sequence(virus_sequence_path):

    virus_number = 34 #virus id +219s
    sequence_all = []
    hash_sequence = []
    for i in range(1,virus_number+1):
        sequence_file = os.path.join(virus_sequence_path,str(i)+".fasta")
        f = open(sequence_file,'r')
        flag = 0
        single_Sequence = []
        for line in f.readlines():
            part_Sequence = line.strip().split("\n")
            part_Sequence_c = [str(i) for i in part_Sequence[0]]
            if flag == 0:
                flag = 1
                continue
            if len(part_Sequence_c) > 0 and part_Sequence_c[0] in ['A','T','G','C']:
                single_Sequence.extend(part_Sequence_c)
        print("ss")
        print(single_Sequence)
        sequence_all.append(single_Sequence)
        hash_sequence.append(hash_seq(single_Sequence,512))
    print(len(sequence_all))
    sequence_all = np.array(sequence_all)
    print(sequence_all.shape)
    #print(hash_sequence)

    new_hash = []
    for i in hash_sequence:
        if len(i) == 512:
            print("yes")
            new_hash.append(i)
        else:
            if len(i) < 512:
                padding = np.zeros(512 - len(i))
                new_hashes = np.concatenate((padding, i), axis=0)
                new_hash.append(new_hashes)
            else:
                new_hash.append(i[0:512])

    new_hash = np.array(new_hash)
    print(np.array(new_hash).shape)
    print(new_hash.reshape(len(hash_sequence), 512))

    np.save("data/HDVD_sequence_hash", new_hash)
    np.save("data/HDVD_id_virus_sequenc", sequence_all)

    return new_hash;






drug_smile_path = os.path.join(HDVD_path, drug_smile)
virus_sequence_path = os.path.join("raw_data/HDVD/","virus_sequence/")
dv_path = os.path.join(HDVD_path, virusdrug)
drug_hash = drug_SMILE_read(drug_smile_path)

#load_interactions(dv_path)

virus_hash = process_virus_sequence(virus_sequence_path)
print(drug_hash.shape)
print(virus_hash.shape)

entity_hash = np.concatenate([drug_hash,virus_hash],axis=0)
print(entity_hash.shape)
np.save("data/entity_hash", entity_hash)