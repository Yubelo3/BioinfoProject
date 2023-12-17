import re
from subprocess import check_output

def tm_score(reference:str,target:str):
    '''
    运行TMalign并从TMalign程序中提取TM-score输出
    '''
    ret = check_output(['./TMalign', reference, target])
    tm_scores = []
    for line in str(ret).split('\\n'):
        if re.match(r'^TM-score=', line):
            score = line.split()[1:2] # Extract the value
            tm_scores.extend(score) # Saving only values
    return float(tm_scores[0])

if __name__=="__main__":
    print(tm_score("pdb/01/201L-A.pdb","pdb/1B/21BI-A.pdb"))