# -*- coding: utf-8 -*-
# @Time    : 2021/11/13 19:13
# @Author  : Bert
# @Email   : 409356793@qq.com
# @File    : WetLabModualation.py
# @Software: PyCharm
# @Version :
# @Content :
import struct
import numpy as np
import random as rd
import pickle
import math
class ConfigDict:
    encodedict = {'00': 'A', '01': 'T', '10': 'C', '11': 'G'}
    revencodedict = {'A': '00', 'T': '01', 'C': '10', 'G': '11'}
    tempUnitLen = 8
    tempUnitStr = '11100100'
    modulecodelenth = 112  ### default length of carrier sequence
    __modulecodestr = ""  # store carrier

    @staticmethod
    def getInforByteNumByModuleCode():
        return (ConfigDict.modulecodelenth)//8

    ##set carrier peroid(unit）
    @staticmethod
    def setModulateUnitCode(moduUnitstr):
        ConfigDict.tempUnitLen = len(moduUnitstr)  # carrier period
        ConfigDict.tempUnitStr = moduUnitstr  # carrier unit
        ConfigDict.__modulecodestr = ""

    #### get carrier sequence
    @staticmethod
    def getModuleCode():
        if ConfigDict.__modulecodestr != "":
            return ConfigDict.__modulecodestr
        kt = (ConfigDict.modulecodelenth)// ConfigDict.tempUnitLen
        tmpstrb = ConfigDict.tempUnitStr * kt
        kt = (ConfigDict.modulecodelenth) % ConfigDict.tempUnitLen
        if kt != 0:
            tmpstrb += ConfigDict.tempUnitStr[:kt]
        return tmpstrb

    ##set carrier
    @staticmethod
    def setModuleCode(lstrb):
        ConfigDict.modulecodelenth = len(lstrb)
        ConfigDict.__modulecodestr = lstrb
        return True

    ### Set the carrier length
    @staticmethod
    def setModuleCodelenth(lenvar):
        ConfigDict.modulecodelenth = lenvar
        return True

    ###Retrieve the carrier length
    @staticmethod
    def getModuleCodelength():
        return ConfigDict.modulecodelenth

##Convert number into binary format with length of bitnum
def ConvertBitSeq(number,bitnum):
    line=bin(number)[2:]
    if len(line)<bitnum:
        line='0'*(bitnum-len(line))+line
    return line

##Convert number into a byte
def ConvertIntByte(number):
    line=bin(number)[2:]
    if len(line)<8:
        line='0'*(8-len(line))+line
    return line

def ConverByteToBaseStr(linestrb,templinestr,size=8):
    line=""
    for i in range(size):
        tmpstr=templinestr[i]
        tmpstr+=linestrb[i]
        line+=ConfigDict.encodedict[tmpstr]
    # line+='\n'
    return line

def ExamineDecodeStr(templineReal,linestrb,templinestr,size=8):
    if templineReal[:size]==templinestr[:size]:
        return linestrb,templineReal
    ###根据模板纠正
    return linestrb,templineReal

def ConverBaseStrToByteList(linebasestr,templinestr):
    tmplist=[]
    templineReal=""
    linestrb=""
    for vt in linebasestr:
        tmp=ConfigDict.revencodedict[vt]
        templineReal+=tmp[0]
        linestrb+=tmp[1]
    linestrb,templineReal=ExamineDecodeStr(templineReal, linestrb, templinestr, len(linestrb))

    for index in range(0,len(linestrb),8):
        byteint=int(linestrb[index:index+8],2)
        tmplist.append(byteint)
    return tmplist

def toBinaryStr(bytevalue,returnstrlen=8):
    tmpstr=str(bin(bytevalue))[2:]
    k=returnstrlen-len(tmpstr)
    final='0'*k
    final+=tmpstr
    return final
def returnStr(linklist):
    line=""
    for vt in linklist:
        vtchr=chr(vt)
        if vtchr=='\n':
            line+='@'
        else:
            line+=vtchr
    return line

########Encode ifilename into ofilename in DNA sequences format
def encodeFile_excluIndex_N(ifilename,ofilename):
    linetemplatestr=ConfigDict.getModuleCode()
    with open(ifilename,"rb") as ifile,open(ofilename,'w',encoding='utf-8') as ofile:
        data=ifile.read()
        ifilesize=len(data)
        lineindex = 0
        for index in range(0,ifilesize,ConfigDict.getInforByteNumByModuleCode()):
            if index+ConfigDict.getInforByteNumByModuleCode()<ifilesize:##可以顺利读取一行
                lineindex+=1
                linearray=data[index:index+ConfigDict.getInforByteNumByModuleCode()]
                linestrb=""
                for ele in linearray:
                    linestrb+=toBinaryStr(ele)
                linestr=ConverByteToBaseStr(linestrb,linetemplatestr,len(linetemplatestr))
                linestr += '\n'
                ofile.write(linestr)
            else:
                linearray = []
                lineindex += 1
                linearray += data[index:]
                linestrb = ""
                for ele in linearray:
                    linestrb += toBinaryStr(ele)
                linestr = ConverByteToBaseStr(linestrb, linetemplatestr, len(linestrb))
                linestr+='\n'
                ofile.write(linestr)

def theta(a, b):
    if a == '-' or b == '-' or a != b:   # gap or mismatch
        return -1
    elif a == b:                         # match
        return 1

def make_score_matrix(seq1, seq2):
    """
    return score matrix and map(each score from which direction)
    0: diagnosis
    1: up
    2: left
    """
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    score_mat = {}
    trace_mat = {}

    for i,p in enumerate(seq1):
        score_mat[i] = {}
        trace_mat[i] = {}
        for j,q in enumerate(seq2):
            if i == 0:                    # first row, gap in seq1
                score_mat[i][j] = -j
                trace_mat[i][j] = 1
                continue
            if j == 0:                    # first column, gap in seq2
                score_mat[i][j] = -i
                trace_mat[i][j] = 2
                continue
            ul = score_mat[i-1][j-1] + theta(p, q)     # from up-left, mark 0
            l  = score_mat[i][j-1]   + theta('-', q)   # from left, mark 1, gap in seq1
            u  = score_mat[i-1][j]   + theta(p, '-')   # from up, mark 2, gap in seq2
            picked = max([ul,l,u])
            score_mat[i][j] = picked
            trace_mat[i][j] = [ul, l, u].index(picked)   # record which direction
    return score_mat, trace_mat

def traceback(seq1, seq2, trace_mat):
    '''
    find one optimal traceback path from trace matrix, return path code
    -!- CAUTIOUS: if multiple equally possible path exits, only return one of them -!-
    '''
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    path_code = ''
    while i > 0 or j > 0:
        direction = trace_mat[i][j]
        if direction == 0:                    # from up-left direction
            i = i-1
            j = j-1
            path_code = '0' + path_code
        elif direction == 1:                  # from left
            j = j-1
            path_code = '1' + path_code
        elif direction == 2:                  # from up
            i = i-1
            path_code = '2' + path_code
    return path_code

def pretty_print_align(seq1, seq2, path_code):
    '''
    return pair alignment result string from
    path code: 0 for match, 1 for gap in seq1, 2 for gap in seq2
    '''
    align1 = ''
    middle = ''
    align2 = ''
    for p in path_code:
        if p == '0':
            align1 = align1 + seq1[0]
            align2 = align2 + seq2[0]
            if seq1[0] == seq2[0]:
                middle = middle + '|'
            else:
                middle = middle + ' '
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        elif p == '1':
            align1 = align1 + '-'
            align2 = align2 + seq2[0]
            middle = middle + ' '
            seq2 = seq2[1:]
        elif p == '2':
            align1 = align1 + seq1[0]
            align2 = align2 + '-'
            middle = middle + ' '
            seq1 = seq1[1:]

    return (align1,align2)

def AlignTwoSeqs(seq1, seq2):
    score_mat, trace_mat = make_score_matrix(seq1, seq2)
    path_code = traceback(seq1, seq2, trace_mat)
    NewSeed,NewCandidate=pretty_print_align(seq1, seq2, path_code)
    return (NewSeed,NewCandidate)

def ConverBaseStrToBinaStr(linebasestr):
    tmplist=[]
    templineReal=""
    linestrb=""
    for vt in linebasestr:
        tmp=ConfigDict.revencodedict[vt]
        templineReal+=tmp[0]
        linestrb+=tmp[1]
    return linestrb,templineReal

def  GroupReads(iflilename,multiplex):
    with open(iflilename,'r',encoding='utf-8') as filept:
        data=filept.readlines()
        data=[vt[:-1] for vt in data]
        encode_num=len(data)//multiplex
        if len(data) %multiplex!=0:
            print("simulation src file has error!")
        cls=[[] for vt in range(encode_num)]
        # for eleid,stline in enumerate(data):
        #     data[eleid]=data[eleid][:-1]
        vcount=0
        for eleid in range(encode_num):
            for tid in range(multiplex):
                cls[eleid].append(vcount)
                vcount+=1
        return cls,data

def ins_del_correct_v2(linestr):
    templatestr = ConfigDict.getModuleCode()
    linestrb,templineRealstrb=ConverBaseStrToBinaStr(linestr)
    errorpos=[]
    seq1, seq2 = AlignTwoSeqs(templineRealstrb, templatestr)
    resline=""
    idvalue=0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != '-':
            resline+=linestr[idvalue]
            errorpos.append('n') # 'n' normal
            idvalue+=1
        elif seq1[i] == seq2[i] and seq1[i] == '-':
            continue
        elif seq1[i] != seq2[i] and seq1[i] == '-':# indicate deletion
             errorpos.append('d')  # 'd' denotes deletion
             if seq2[i]=='0':
                 if np.random.rand()<=0.5:
                     resline+='A'
                 else:
                     resline+='T'
             else:
                 if np.random.rand()<=0.5:
                     resline+='G'
                 else:
                     resline+='C'
        elif seq1[i] != seq2[i] and seq2[i] == '-':## insertion
              idvalue+=1
        elif seq1[i] != seq2[i]: # substitute
            errorpos.append('s')  #
            if linestr[idvalue] in ['A','T']:
                if np.random.rand() <= 0.5:
                    resline += 'G'
                else:
                    resline += 'C'
            else:
                if np.random.rand() <= 0.5:
                    resline += 'A'
                else:
                    resline += 'T'
            idvalue+=1
    return resline,errorpos
def generateCandidate_1(newcontent,errorarray):
    tsize=len(newcontent)
    presline=""
    for i in range(len(newcontent[0])):
        aplabet={'A':0,'T':0,'G':0,'C':0}
        linsr=''
        for vx in range(tsize):
            if errorarray[vx][i]=='n':
                aplabet[newcontent[vx][i]]+=1
                linsr+=newcontent[vx][i]
        presline+=sorted(aplabet.items(), key = lambda kv:(kv[1], kv[0]),reverse=True)[0][0]
    return presline

#### error correction of reads for each cluster and output consensus sequences
def cluster_indel_correction(clusterdict,data):
    rectiyreads=[]
    errorstatarray=[]
    ####indel error correction...
    for dtline in data:
        tuele=ins_del_correct_v2(dtline)
        rectiyreads.append(tuele[0])
        errorstatarray.append(tuele[1])
    # rectiyreads=data
    ####consensus sequence...
    content = []
    for clusele in clusterdict:
        if len(clusele)!=0:
            tmpcontent = [rectiyreads[ikd] for ikd in clusele]
            tmperray = [errorstatarray[ikd] for ikd in clusele]
            cite = generateCandidate_1(tmpcontent, tmperray)
            content.append(cite)
    return content

####cluster and error correction
def cluster_majority(ifilename,ofilename,multiplex):
    cls,data=GroupReads(ifilename,multiplex) ## In practical, this function should have clustering function such as Starcode...
    ###Error correction for each cluster##########
    rectifyreads=cluster_indel_correction(cls,data)
    with open(ofilename,"w",encoding='utf-8')as wfile:
        for line in rectifyreads:
            wfile.write(line)
            wfile.write("\n")

def decodefile_excluIndex(ifilename,ofilename):
    linetemplatestr = ConfigDict.getModuleCode()
    with open(ifilename, "r",encoding='utf-8') as ifile, open(ofilename, 'wb') as ofile:
        data=ifile.readlines()
        for linebasestr in data:
            tmplist=ConverBaseStrToByteList(linebasestr[:-1], linetemplatestr)
            for vt in tmplist:
                ain = struct.pack('B', vt)
                ofile.write(ain)
    # if os.path.exists(ifilename):
    #     os.remove(ifilename)

#random return a base different from ch
def returnBase(ch):
    t=[chv for chv in ['A','T','G','C'] if chv!=ch]
    return t[rd.randint(0,2)]

## statistical character frequency of the orginal encoding text,
# this output of this function is used to compare the data recovery
def statisTemplateFile(filename):
    chrdict=dict()
    with open(filename,"rb") as file,open("template_refer_lib.lib","wb") as ft:
        data = file.read()
        ifilesize = len(data)
        contentfile = ""
        for ik in range(len(data)):
            contentfile += chr(data[ik])
        # content=file.readlines()
        # content=''.join(content)
        content=contentfile
        charset=set(content)
        chrset=list(charset)
        for chele in chrset:
            chrdict[chele]=0
        for chele in content:
            chrdict[chele]+=1
        pickle.dump(chrdict,ft)
        pickle.dump(len(content),ft)

##### Generate Simulation sequenced reads according to the encoded sequence file (ifilename)
###### wtih specified error rate , sequence copies (multipleX), and repeat_id（ Repeatid).
####### The error rate is the sum of insertion rate, deletion rate, and substitution rate,and all these error rate is equal.
def  StimulateBaseError_equalCopies(ifilename,errorrate,multipleX,Repeatid):
    ofilename=ifilename.replace(".txt","")
    ofilename+="_X"+str(multipleX)
    ofilename+="_R"+str(Repeatid)
    ofilename+="_E"+str(errorrate)
    oifilename=ofilename #存放行实际对应的行
    ofilename+=".txt"
    oifilename+="rowindex.txt"
    # matofilename = 'template_'+oifilename
    # matofilename = matofilename.replace('.txt', '.mat')
    matrix=[]
    with open(ifilename,"r",encoding="UTF-8") as file:
        arrayline=file.readlines()
    lineIndexList=[]
    for kindex in range(len(arrayline)):
        for multiid in range(multipleX):
            lineIndexList.append(kindex)
    errorMatrxi=[]
    for kindex in lineIndexList:
        # kindex=rd.randint(0,len(arrayline)-1)
        matrix.append([list(arrayline[kindex]),[kindex]])

    for index in range(len(matrix)):
        errorlist=np.random.rand((len(matrix[index][0])-1))
        substitulist=[]
        inslist=[]
        dellist=[]
        for kid in range(len(errorlist)):
            if errorlist[kid]<=errorrate:
                krd = np.random.randint(0, 3)
                if krd == 0:
                    substitulist.append(kid)
                elif krd==1:
                    inslist.append(kid)
                else:
                    dellist.append(kid)
        errorMatrxi.append([substitulist,inslist,dellist])
    #     set substitution error
        for tx in substitulist:
            matrix[index][0][tx] = returnBase(matrix[index][0][tx])
            matrix[index][1].append((tx,'S'))
     #   set deletion error
        for tx in dellist:
            matrix[index][0][tx] = 'D'
            matrix[index][1].append((tx,'D'))
        # set insertion error
        inslist=sorted(inslist,reverse=True)
        basechrset = ['A', 'G', 'T', 'C']
        for tx in inslist:
            matrix[index][0].insert(tx,basechrset[np.random.randint(0,4)])
            matrix[index][1].append((tx, 'I'))

    ordlist=[x for x in range(len(matrix))]

    with open(ofilename,"w",encoding="UTF-8") as file1,open(oifilename,"w",encoding="UTF-8") as file2:
       # inferInforArray=[]
       for pid in ordlist:
           file1.write(''.join(matrix[pid][0]).replace('D',''))
           speflag=True
           file2.write("{} ".format(matrix[pid][1][0]))
           tmpvarlist=matrix[pid][1][1:]
           tmlist=sorted(tmpvarlist,key=lambda xt:xt[0])
           for eid in tmlist:
               file2.write("({},{}) , ".format(eid[0], eid[1]))
           file2.write("\n")
           # inferInforArray.append(tmlist)
    return (ofilename, oifilename)

def isValidErrorPos(ikt,indexarray,consecIndelNum):
    flagvar=True
    if ikt+consecIndelNum-1>=len(indexarray):
        return False
    for ikx in range(consecIndelNum):
        if indexarray[ikt+ikx]!=0:
            flagvar=False
            break
    return flagvar

def comFileLikely1(filename1):### Compare the similarity between decoded file and original file
    chrdict=dict()
    filesize=0
    with open("template_refer_lib.lib","rb") as file1:
        chrdict=pickle.load(file1)
        filesize=pickle.load(file1)
    with open(filename1, "rb") as file:
        data = file.read()
        contentfile=""
        for ik in range(len(data)):
            contentfile+=chr(data[ik])
    count=0
    for vch in contentfile:
        if vch in chrdict and chrdict[vch]>0:
            chrdict[vch]-=1
        else:#decode a character not present in original text
            count+=1
    sum=0
    for vch in chrdict:
        if chrdict[vch]>0:
            sum+=abs(chrdict[vch])
    # sum+=count
    return 1-sum/filesize


def TestDemo(modulateCodeLength,error,seqdepth):
    ConfigDict.setModuleCodelenth(modulateCodeLength)
    ink=1 ## carrier id
    basefile='encodebase_{}_.txt'.format(ink)
    print("The {} carrier".format(ink))

    #### encode the the file 'grandmother.txt‘
    encodeFile_excluIndex_N('grandmother.txt', basefile)

    ######## Statistic the frequency of characters in the 'grandmother.txt', this function is used in the final accuracy of recovery data.
    statisTemplateFile('grandmother.txt')

    ###According to the basefile, simultate the synthesis and sequencing process.
    StimulateBaseError_equalCopies(basefile, error, seqdepth, 0)
    ofilename = basefile.replace(".txt", "")
    ofilename += "_X" + str(seqdepth)
    ofilename += "_R" + str(0)
    ofilename += "_E" + str(error)
    ofilename += ".txt"

    #### modulation-based error detection and error correction
    cluster_majority(ofilename, ofilename.replace('.txt', '_consensus.txt'),seqdepth)

    ####decoding the consensus sequences according to the modulation rule.
    decodefile_excluIndex(ofilename.replace('.txt', '_consensus.txt'),
                          ofilename.replace('.txt', '_result.txt'))
    k = comFileLikely1(ofilename.replace('.txt', '_result.txt'))

    print('modulatecodeid,{},errorate,{},seqdepth,{},Data recovery accuracy,{}'.format(ink,error, seqdepth, k))

if __name__=="__main__":

    ######################## Simulation Error correction performance.##########################
    ## Notes:  carrier signal consists of substirng '11100100',and the length of the encoded sequence is 112.
    TestDemo(112,0.4,80)

    # TestDemo(112, 0.4, 60)
