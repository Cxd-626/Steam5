import itertools
import argparse
import os,sys,re
import numpy as np 
from Bio import SeqIO
import pickle as pkl
import scipy.io as sio

def CalculateMatrix(data, order, k):
    if k == 1:
        matrix = np.zeros((len(data[0]), 4))
        for i in range(len(data[0]) - 1): # position
            for j in range(len(data)):
                matrix[i][order[data[j][i:i+1]]] += 1
    elif k == 2:
        matrix = np.zeros((len(data[0]) - 1, 16))
        for i in range(len(data[0]) - 2): # position
            for j in range(len(data)):
                matrix[i][order[data[j][i:i+2]]] += 1
    else:
        matrix = np.zeros((len(data[0]) - 2, 64))
        for i in range(len(data[0]) - 2): # position
            for j in range(len(data)):
                matrix[i][order[data[j][i:i+3]]] += 1           
    return matrix

def test_PSP(train_datapath, test_datapath, k):  
    
    sequences  = read_fasta(train_datapath)
    train_positive = []
    for i in range(int(len(sequences)/2)):
        train_positive.append(str(sequences[i]))
        
    train_negative = []
    for j in range(int(len(sequences)/2),int(len(sequences))):
        train_negative.append(str(sequences[j]))
        
    train_p_num = len(train_positive)
    train_n_num = len(train_negative)
    
    test_seq = read_fasta(test_datapath)
    test = []
    
    for pos in test_seq:
        test.append(str(pos))
        
    test_num = len(test)
    test_l = len(test[0])

    nucleotides = ['A', 'C', 'G', 'U']
    
    if k == 1 :
        nuc = [n1 for n1 in nucleotides]
        order = {}
        for i in range(len(nuc)):
            order[nuc[i]] = i
        
        matrix_po = CalculateMatrix(train_positive, order, 1)
        matrix_ne = CalculateMatrix(train_negative, order, 1)

        F1 = matrix_po/train_p_num
        F2 = matrix_ne/train_n_num       
        F = F1 - F2

        code = []
        for sequence in test:  
            for j in range(len(sequence)):                
                number = F[j][order[sequence[j:j+1]]]
                code.append(number)  
        code = np.array(code)
        code = code.reshape((test_num,test_l))
        
    
    elif k == 2:
        dnuc = [n1 + n2  for n1 in nucleotides for n2 in nucleotides]
        order = {}
        for i in range(len(dnuc)):
            order[dnuc[i]] = i
        
        matrix_po = CalculateMatrix(train_positive, order, 2)
        matrix_ne = CalculateMatrix(train_negative, order, 2)
        
        F1 = matrix_po/train_p_num
        F2 = matrix_ne/train_n_num       
        
        F = F1 - F2
        
        code = []
        for sequence in test:  
            for j in range(len(sequence)-1):                
                number = F[j][order[sequence[j:j+2]]]
                code.append(number)  
        code = np.array(code)
        code = code.reshape((test_num,test_l-1))
        
       
        
    else:
        tnuc = [n1 + n2 + n3 for n1 in nucleotides for n2 in nucleotides for n3 in nucleotides]
        order = {}
        for i in range(len(tnuc)):
            order[tnuc[i]] = i
        
        matrix_po = CalculateMatrix(train_positive, order, 3)
        matrix_ne = CalculateMatrix(train_negative, order, 3)
        
        F1 = matrix_po/train_p_num
        F2 = matrix_ne/train_n_num       
        
        F = F1 - F2
    
        code = []
        for sequence in test:  
            for j in range(len(sequence)-2):                
                number = F[j][order[sequence[j:j+3]]]
                code.append(number)  
        code = np.array(code)
        code = code.reshape((test_num,test_l-2))
        
    return code



def Kmer(sequences):
    AA = 'ACGU'
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    Kmer_feature = []

    for seq in sequences:
        kmer1 = [0] * 4
        for j in range(len(seq)):
            kmer1[AADict[seq[j]]] = kmer1[AADict[seq[j]]] + 1
        if sum(kmer1) != 0:
            kmer1 = [i / sum(kmer1) for i in kmer1]

        kmer2 = [0] * 16
        for j in range(len(seq) - 2 + 1):
            kmer2[AADict[seq[j]] * 4 + AADict[seq[j + 1]]] = kmer2[AADict[seq[j]] * 4 + AADict[seq[j + 1]]] + 1
        if sum(kmer2) != 0:
            kmer2 = [i / sum(kmer2) for i in kmer2]

        kmer3 = [0] * 64
        for j in range(len(seq) - 3 + 1):
            kmer3[AADict[seq[j]] * 16 + AADict[seq[j + 1]] * 4 + AADict[seq[j + 2]]] = kmer3[AADict[seq[j]] * 16 + AADict[seq[j + 1]] * 4 + AADict[seq[j + 2]]] + 1
        if sum(kmer3) != 0:
            kmer3 = [i / sum(kmer3) for i in kmer3]

        kmer = kmer1 + kmer2 + kmer3
        Kmer_feature.append(kmer)
    return Kmer_feature 
 
myDiIndex = {
    'AA': 0, 'AC': 1, 'AG': 2, 'AU': 3,
    'CA': 4, 'CC': 5, 'CG': 6, 'CU': 7,
    'GA': 8, 'GC': 9, 'GG': 10, 'GU': 11,
    'UA': 12, 'UC': 13, 'UG': 14, 'UU': 15
}

baseSymbol = 'ACGU'


def get_kmer_frequency(sequence, kmer):
    myFrequency = {}
    for pep in [''.join(i) for i in list(itertools.product(baseSymbol, repeat=kmer))]:
        myFrequency[pep] = 0
    for i in range(len(sequence) - kmer + 1):
        myFrequency[sequence[i: i + kmer]] = myFrequency[sequence[i: i + kmer]] + 1
    for key in myFrequency:
        myFrequency[key] = myFrequency[key] / (len(sequence) - kmer + 1)
    return myFrequency


def correlationFunction(pepA, pepB, myIndex, myPropertyName, myPropertyValue):
    CC = 0
    for p in myPropertyName:
        CC = CC + (float(myPropertyValue[p][myIndex[pepA]]) - float(myPropertyValue[p][myIndex[pepB]])) ** 2
    return CC / len(myPropertyName)


def correlationFunction_type2(pepA, pepB, myIndex, myPropertyName, myPropertyValue):
    CC = 0
    for p in myPropertyName:
        CC = CC + float(myPropertyValue[p][myIndex[pepA]]) * float(myPropertyValue[p][myIndex[pepB]])
    return CC


def get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, kmer):
    thetaArray = []
    for tmpLamada in range(lamadaValue):
        theta = 0
        for i in range(len(sequence) - tmpLamada - kmer):
            theta = theta + correlationFunction(sequence[i:i + kmer],
                                                sequence[i + tmpLamada + 1: i + tmpLamada + 1 + kmer], myIndex,
                                                myPropertyName, myPropertyValue)
        thetaArray.append(theta / (len(sequence) - tmpLamada - kmer))
    return thetaArray


def get_theta_array_type2(myIndex, myPropertyName, myPropertyValue, lamadaValue, sequence, kmer):
    thetaArray = []
    for tmpLamada in range(lamadaValue):
        for p in myPropertyName:
            theta = 0
            for i in range(len(sequence) - tmpLamada - kmer):
                theta = theta + correlationFunction_type2(sequence[i:i + kmer],
                                                          sequence[i + tmpLamada + 1: i + tmpLamada + 1 + kmer],
                                                          myIndex,
                                                          [p], myPropertyValue)
            thetaArray.append(theta / (len(sequence) - tmpLamada - kmer))
    return thetaArray


def PCPseDNC(sequences):
    
    myPropertyName = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'A-philicity', 'Propeller twist',
                'Duplex stability:(freeenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
                'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH', 'Electron_interaction',
                'Hartman_trans_free_energy', 'Helix-Coil_transition', 'Lisser_BZ_transition', 'Polar_interaction',
                'SantaLucia_dG', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Sugimoto_dG', 'Sugimoto_dH',
                'Sugimoto_dS', 'Duplex tability(disruptenergy)', 'Stabilising energy of Z-DNA', 'Breslauer_dS',
                'Ivanov_BA_transition', 'SantaLucia_dH', 'Stacking_energy', 'Watson-Crick_interaction',
                'Dinucleotide GC Content', 'Twist', 'Tilt', 'Roll', 'Shift', 'Slide', 'Rise']
    dataFile = 'Phychepro.data'
    with open(dataFile,'rb') as f:
        myPropertyValue = pkl.load(f)
    lamadaValue = 2 #20
    weight =  0.1 #0.9
    myIndex = myDiIndex
    PCPseDNC_feature = []
    for i in sequences:
        code = []
        dipeptideFrequency = get_kmer_frequency(i, 2)
        thetaArray = get_theta_array(myIndex, myPropertyName, myPropertyValue, lamadaValue, i, 2)
        for pair in sorted(myIndex.keys()):
            code.append(dipeptideFrequency[pair] / (1 + weight * sum(thetaArray)))
        for k in range(17, 16 + lamadaValue + 1):
            code.append((weight * thetaArray[k - 17]) / (1 + weight * sum(thetaArray)))
        PCPseDNC_feature.append(code)
    return PCPseDNC_feature
 
#PseEIIP
 
def TriNcleotideComposition(sequence, base):
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    tnc_dict = {}
    for triN in trincleotides:
        tnc_dict[triN] = 0
    for i in range(len(sequence) - 2):
        tnc_dict[sequence[i:i + 3]] += 1
    for key in tnc_dict:
       tnc_dict[key] /= (len(sequence) - 2)
    return tnc_dict

def PseEIIP(fastas):
    base = 'ACGU'
    EIIP_dict = {
        'A': 0.1260,
        'C': 0.1340,
        'G': 0.0806,
        'U': 0.1335 }
    trincleotides = [nn1 + nn2 + nn3 for nn1 in base for nn2 in base for nn3 in base]
    EIIPxyz = {}
    for triN in trincleotides:
        EIIPxyz[triN] = EIIP_dict[triN[0]] + EIIP_dict[triN[1]] + EIIP_dict[triN[2]]
    encodings = []
    for sequence in fastas:
        code = []
        trincleotide_frequency = TriNcleotideComposition(sequence, base)
        code = code + [EIIPxyz[triN] * trincleotide_frequency[triN] for triN in trincleotides]
        encodings.append(code)
    return encodings



def read_fasta(datapath):
    if os.path.exists(datapath) == False:
        print('Error: file " %s " does not exist.' % datapath)
        sys.exit(1)
    with open(datapath) as f:
        record = f.readlines()
    if re.search('>',record[0]) == None:
        print('Error: the input file " %s " must be fasta format!' % datapath)
        sys.exit(1)
        
    sequences = list(SeqIO.parse(datapath, "fasta"))
    sequence = []
    for i in range(len(sequences)):
        sequence.append(str(sequences[i].seq))
    return sequence


def extract_features(datapath,speciesfile):
    if speciesfile == 'Arabidopsis':
        train_datapath = "./data/Arabidopsis/Arabidopsis_train.fasta"
        PSNP = np.array(test_PSP(train_datapath, datapath, 1))
        PSDP = np.array(test_PSP(train_datapath, datapath, 2))
        PSTP = np.array(test_PSP(train_datapath, datapath, 3))
        PSP = np.concatenate((PSNP, PSDP, PSTP), axis=1)
        seq = read_fasta(datapath)
        Kmer1 = np.array(Kmer(seq))
        PCPseDNC1 = np.array(PCPseDNC(seq))
        PseEIIP1 = np.array(PseEIIP(seq))
        feature_vector = np.concatenate((PSP, Kmer1, PCPseDNC1, PseEIIP1), axis=1)
    elif speciesfile == 'Mouse':
        train_datapath = "./data/Mouse/mouse_train.fasta"
        PSNP = np.array(test_PSP(train_datapath, datapath, 1))
        PSDP = np.array(test_PSP(train_datapath, datapath, 2))
        PSTP = np.array(test_PSP(train_datapath, datapath, 3))
        PSP = np.concatenate((PSNP, PSDP, PSTP), axis=1)
        seq = read_fasta(datapath)
        Kmer1 = np.array(Kmer(seq))
        PCPseDNC1 = np.array(PCPseDNC(seq))
        PseEIIP1 = np.array(PseEIIP(seq))
        pri_feature_vector = np.concatenate((PSP, Kmer1, PCPseDNC1, PseEIIP1), axis=1)
        f = sio.loadmat(r'Fscore.mat')
        score  = f['F1']
        feature_vector = []
        for i in score[1,0:185]:
            feature_vector.append(pri_feature_vector[:,int(i-1)])
        feature_vector = np.array(feature_vector).T
    return np.array(feature_vector)
        
    
        
def main():
    parser = argparse.ArgumentParser(description='Staem5: a novel stacked ensemble method for prediction of m5C site')
    parser.add_argument('--input',dest='inputpath',type=str,required=True, help='query RNA positive sequences to be predicted in fasta format.')
    parser.add_argument('--species', dest='speciesfile', type=str, required=False, help='--species indicates the specific species, currently we accept \'Arabidopsis\' or \'Mouse\'', default=None)
    parser.add_argument('--output',dest='outputfile',type=str,required=False, help='save the prediction results in csv format.')
    args = parser.parse_args()
    
    inputpath = args.inputpath
    outputfile = args.outputfile
    speciesfile = args.speciesfile
    
    outputfile_original = outputfile
    
    if outputfile_original == None:
        outputfile_original = ''
    try:
        if outputfile_original == None:
            outputfile = 'output'
        outputfile = outputfile_original + '_' + speciesfile    
        
        if speciesfile == 'Arabidopsis':
            vector  = extract_features(inputpath,speciesfile) 
            model = pkl.load(open("./model/Arabidopsis_stackClf.pickle.model", "rb"))
            predictions = model.predict_proba(vector)
            sequence = read_fasta(inputpath)
            seq = []
            for i in sequence:
                seq.append(str(i))
            probability = ['%.5f' % float(i) for i in predictions[:,1]]
            with open(outputfile,'w') as f:
                for i in range(int(len(vector))):
                    if float(probability[i]) > 0.5:
                        f.write(probability[i]+ '*' + '\t')
                        f.write(seq[i] + '*' + '\t')
                        f.write('1' + '\n')
                    else:
                        f.write(probability[i] + '*' + '\t')
                        f.write(seq[i] + '*' + '\t')
                        f.write('0' + '\n')
            print('output are saved in ' + outputfile + ', and those with probability greater than 0.5 are marked with *')
       
        elif speciesfile == 'Mouse':
            vector  = extract_features(inputpath,speciesfile) 
            model = pkl.load(open("./model/Mouse_stackClf.pickle.model", "rb"))
            predictions = model.predict_proba(vector)
            sequence = read_fasta(inputpath)
            seq = []
            for i in sequence:
                seq.append(str(i))
            probability = ['%.5f' % float(i) for i in predictions[:,1]]
            with open(outputfile,'w') as f:
                for i in range(int(len(vector))):
                    if float(probability[i]) > 0.5:
                        f.write(probability[i]+ '*' + '\t')
                        f.write(seq[i] + '*' + '\t')
                        f.write('1' + '\n')
                    else:
                        f.write(probability[i] + '*' + '\t')
                        f.write(seq[i] + '*' + '\t')
                        f.write('0' + '\n')
            print('output are saved in ' + outputfile + ', and those with probability greater than 0.5 are marked with *')
    
    except Exception as e:
        print('Please check the format of your predicting data!')
        sys.exit(1)
        
if __name__ == "__main__":
    main()
 



