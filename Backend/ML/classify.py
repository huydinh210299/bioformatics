## import thu vien
import os
from Pfeature.pfeature import paac_wp
import pandas as pd
import pickle
from Bio import SeqIO
import random
import string
import threading


class Classify:

  # Ham thuc hien generate du lieu
  def paac(seft, input):
    a = input.rstrip('fasta')
    output = a + 'paac.csv'
    df_out = paac_wp(input, output,4,0.2)
    df = pd.read_csv(output)

  # Ham luu seq truyen vao thanh file va genarate ra du lieu
  def paacBySeq(seft, seq):
    fileName = string.ascii_lowercase
    fileName = ''.join(random.choice(fileName) for i in range(10))
    fileName = fileName + ".fasta"
    seft.filename = fileName
    with open(fileName, "w") as file:
      # Writing data to a file
      file.write(">001\n")
      file.write(seq)
    genarateThread = threading.Thread(target=seft.paac,args=(fileName,))
    genarateThread.start()
    return fileName

  def paacByFile(seft, input_fname):
    fileName = string.ascii_lowercase
    fileName = ''.join(random.choice(fileName) for i in range(10))
    fileName = fileName + ".fasta"
    os.rename(input_fname, fileName)
    genarateThread = threading.Thread(target=seft.paac,args=(fileName,))
    genarateThread.start()
    return fileName

  # Ham load model
  def loadModel(seft):
    with open('./model/model.pkl', 'rb') as file:
      pickle_model = pickle.load(file)
    return pickle_model 

  # Ham chuyen tu file fasta thanh array
  def fasta_to_Array(seft, input):
    arr = []
    for record in SeqIO.parse(input, "fasta"):
      item = [record.id,str(record.seq)]
      arr.append(item)
    return arr

  # Ham thuc hien phan loai voi file
  def antiMicroable(seft, input_fname):
    df = pd.read_csv(input_fname)
    model = seft.loadModel()
    yPred = model.predict(df)
    yPred = pd.DataFrame(yPred,columns=['class'])
    print(input_fname)
    fasta_file = input_fname.split(".")[1]
    fasta_file = '.' + fasta_file + '.fasta'
    print(fasta_file)
    mydict = seft.fasta_to_Array(fasta_file)
    df1 = pd.DataFrame(mydict,columns=['ID', 'Seq'])
    result = pd.concat([df1, yPred], axis=1, ignore_index=True)
    result.columns = ["ID", "Seq", "Class"]
    return result
