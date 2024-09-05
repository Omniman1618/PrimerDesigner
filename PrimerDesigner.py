from __future__ import print_function
from base64 import b64encode
import requests
import json
from urllib import request, parse
import pandas as pd
import winsound


def get_access_token(client_id, client_secret, idt_username, idt_password):
    """
    Create the HTTP request, transmit it, and then parse the response for the 
    access token.
    
    The body_dict will also contain the fields "expires_in" that provides the 
    time window the token is valid for (in seconds) and "token_type".
    """

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = { "Content-Type" : "application/x-www-form-urlencoded",
                        "Authorization" : "Basic " + authorization_string }
                    
    data_dict = {   "grant_type" : "password",
                    "scope" : "test",
                    "username" : idt_username,
                    "password" : idt_password }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token", 
                                    data = request_data, 
                                    headers = request_headers,
                                    method = "POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()
    
    # Error and return the response from the endpoint if there was a problem
    if (response.status != 200):
        raise RuntimeError("Request failed with error code:" + response.status + "\nBody:\n" + body)
    
    body_dict = json.loads(body)
    return body_dict["access_token"]
        
#Define function used to generate complement
def flip_seq(init_seq):
    flipped_seq = ""
    seq_len = len(init_seq)
    while seq_len > 0:
        stripped_base = init_seq[-1]
        flipped_seq = flipped_seq + stripped_base
        init_seq = init_seq[:-1]
        seq_len = seq_len -1
    
    flipped_seq = flipped_seq.replace("C","g")
    flipped_seq = flipped_seq.replace("G","c")
    flipped_seq = flipped_seq.replace("A","t")
    flipped_seq = flipped_seq.replace("T","a")
    flipped_seq = flipped_seq.upper()


    return flipped_seq    

def PrimerGen(BindingArea):
    PrimerOptions = []
    ba_len = len(BindingArea)
    prim_len_start = min_primer_len

    while prim_len_start < max_primer_len + 1:
        num_of_options = ba_len - prim_len_start + 1
        ps = 0
        pf = prim_len_start
        pl = pf - ps
        
        while num_of_options > 0:
            primer_option = BindingArea[ps:pf]
            gc_count = primer_option.count("G") + primer_option.count("C")
            primer_gc_percent = gc_count / pl
            occurences = whole_sequence.count(primer_option) + flipped_whole_seq.count(primer_option)
            if primer_gc_percent >= min_gc_percent and occurences == 1:
                PrimerOptions.append(primer_option)

            else:
                None

            ps = ps + 1
            pf = pf + 1
            num_of_options = num_of_options -1

        prim_len_start = prim_len_start + 1

    return PrimerOptions

def homodimer(bearer, sequence):
    url = f'https://www.idtdna.com/restapi/v1/OligoAnalyzer/HeteroDimer?primary={sequence}&secondary={sequence}'

    payload = {}
    headers = {
    'Authorization': f'Bearer {bearer}',
    'Cookie': 'ARRWestffinity=6b8cb0020cc3e79eb7829ae732e01f82a8e045721802c02a438cbe1d24185cb6'
    }

    response = requests.request("POST", url, headers=headers, data=payload)
    response = response.json()
    delta_g = response[0]["DeltaG"]
    return delta_g

def heterodimer(bearer, forward_sequence, reverse_sequence):
    url = f'https://www.idtdna.com/restapi/v1/OligoAnalyzer/HeteroDimer?primary={forward_sequence}&secondary={reverse_sequence}'

    payload = {}
    headers = {
    'Authorization': f'Bearer {bearer}',
    'Cookie': 'ARRWestffinity=6b8cb0020cc3e79eb7829ae732e01f82a8e045721802c02a438cbe1d24185cb6'
    }

    response = requests.request("POST", url, headers=headers, data=payload)
    response = response.json()
    delta_g = response[0]["DeltaG"]
    return delta_g

def analyze_primer(bearer, sequence):
    url = "https://www.idtdna.com/restapi/v1/OligoAnalyzer/Analyze"

    payload = json.dumps({
    "Sequence": sequence,
    "NaConc": 50,
    "MgConc": 3,
    "DNTPsConc": 0.8,
    "OligoConc": 0.2,
    "NucleotideType": "DNA"
    })
    headers = {
    'Content-Type': 'application/json',
    'Authorization': f'Bearer {bearer}',
    'Cookie': 'ARRWestffinity=6b8cb0020cc3e79eb7829ae732e01f82a8e045721802c02a438cbe1d24185cb6'
    }

    response = requests.request("POST", url, headers=headers, data=payload)
    return response
   
def primer_analyze_homo(bearer, sequence):
    deltaG = homodimer(bearer, sequence)
    response = analyze_primer(bearer, sequence)
    response = response.json()
    analyzer_values = response

    sequence = sequence
    TM = analyzer_values["MeltTemp"]
    length = analyzer_values["Length"]
    

    """
    Paramters Used

    min_primer_len = 15
    max_primer_len = 25
    min_dimer_value = -10
    TM_min = 58
    TM_max = 62
    min_gc_percent = 0.50
    """
    
    if (deltaG > min_dimer_value) and (TM_max > TM > TM_min):

        candidate = pd.DataFrame({
            "sequence": [sequence],
            "TM": [TM],
            "homodimer": [deltaG],
            "length": [length]
        })
        
        
        print(candidate)
        return candidate

        
    else:
        return "Skip"
        
#Gets bearer token   
if __name__ == "__main__":
    client_id = "INSERT INFO HERE"
    client_secret = "INSERT INFO HERE"
    idt_username = "INSERT INFO HERE"
    idt_password = "INSERT INFO HERE"
    
    token = get_access_token(client_id, client_secret, idt_username, idt_password)
    #print(token)

#User inputs genetic code~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Generates all primer options within length requirements and GC content requirments
#sliced string = string[start:end:step]


whole_sequence =input("Paste entire sequence in here:\n").upper()
fpba = input("Paste length of DNA where forward primer will bind:\n").upper()
rpba_raw = input("Paste length of DNA where reverse primer will bind:\n").upper()

#Default parameters
min_primer_len = 19
max_primer_len = 21
min_dimer_value = -10
TM_min = 58
TM_max = 62
min_gc_percent = 0.50

#Display Default Parameters
print("Default Settings: ")
print("min_primer_len: ", min_primer_len)
print("max_primer_len: ", max_primer_len)
print("min_dimer_value: ", min_dimer_value)
print("TM_min: ", TM_min)
print("TM_max: ", TM_max)
print("min_gc_percent: ", min_gc_percent)

#Have user either accept or insert new default parameters
keep_def_settings = input("Keep Default Settings? y/n:\n")

if keep_def_settings == "y":
    print("Default Settings Kept")
else:
    min_primer_len = int(input("min_primer_len:\n"))
    max_primer_len = int(input("max_primer_len:\n"))
    min_dimer_value = float(input("min_dimer_value:\n"))
    TM_min = float(input("TM_min:\n"))
    TM_max = float(input("TM_max:\n"))
    min_gc_percent = float(input("min_gc_percent:\n"))

#Flips whole sequence and the reverse primer binding area to the 3' ot 5' orientation
flipped_whole_seq = flip_seq(whole_sequence)
rpba = flip_seq(rpba_raw)

#Generates and prints all candidates within the acceptable GC content and length parameters
fp_init_candidates = PrimerGen(fpba)
rp_init_candidates = PrimerGen(rpba)

print(fp_init_candidates)
print(rp_init_candidates)

#Defines the empty dataframes for the candidates
fp_analyzed_candidates = pd.DataFrame({
    "sequence": [],
    "TM": [],
    "homodimer": [],
    "length": []
})

rp_analyzed_candidates = pd.DataFrame({
    "sequence": [],
    "TM": [],
    "homodimer": [],
    "length": []
})

#Analyzes/checks the homodimer values of all of the initial candidates. Eliminates those that do not meet the parameters, and then sorts based off of homodimer. Prints the sorted lists
print("Analyzing forward primers...........................................")
for fp_init_candidate in fp_init_candidates:
    candidate = primer_analyze_homo(token, fp_init_candidate)
    if len(candidate) > 0 and isinstance(candidate, pd.DataFrame):
        fp_analyzed_candidates = pd.concat([fp_analyzed_candidates, candidate], ignore_index=True)

sorted_fp_candidates = fp_analyzed_candidates.sort_values(by='homodimer', ascending=False)
print("Analyzing reverse primers...........................................")
for rp_init_candidate in rp_init_candidates:
    candidate = primer_analyze_homo(token, rp_init_candidate)
    if len(candidate) > 0 and isinstance(candidate, pd.DataFrame):
        rp_analyzed_candidates = pd.concat([rp_analyzed_candidates, candidate], ignore_index=True)

sorted_rp_candidates = rp_analyzed_candidates.sort_values(by='homodimer', ascending=False)


print("fp candidates: \n",sorted_fp_candidates)
print("rp candidates: \n", sorted_rp_candidates)
print(sorted_rp_candidates.iloc[0]["sequence"])

#Counts the number of forward and reverse primer candidates, then sets the default variables for heterodimer analysis
fp_candidate_count = sorted_fp_candidates.shape[0]
rp_candidate_count = sorted_rp_candidates.shape[0]


rp_number = 0
fp_number = 0
next_switch = "fp"
primary = "fp"
one_done = 0

#Loops through the primers, alternating between FP and RP for primary primer to find the best pair with an acceptable heterodimer
while True:
    if primary == "fp":
            cycles_to_run = fp_number

    elif primary == "rp":
        cycles_to_run = rp_number
         
    #Loop through one primary and all secondaries that are equal or better in rank
    loop_num = 0
    #cycles_to_run starts out at zero. When it finishes the last one and loops back, it will kill the loop and move on to the next one
    while loop_num <= cycles_to_run:
        if primary == "fp":
            internal_rp_num = loop_num
            internal_fp_num = fp_number
            fp_hetero_anal = sorted_fp_candidates.iloc[internal_fp_num]["sequence"]
            rp_hetero_anal = sorted_rp_candidates.iloc[internal_rp_num]["sequence"]
        else: 
            internal_rp_num = rp_number
            internal_fp_num = loop_num
            fp_hetero_anal = sorted_fp_candidates.iloc[internal_fp_num]["sequence"]
            rp_hetero_anal = sorted_rp_candidates.iloc[internal_rp_num]["sequence"]
        
        het_deltaG = heterodimer(token, fp_hetero_anal, rp_hetero_anal)

        if het_deltaG > min_dimer_value:
            break

        else:
            loop_num += 1
    else:
        if one_done == 0:
            if primary == "fp":
                primary = "rp"
                fp_number += 1
                if fp_number == fp_candidate_count:
                    one_done += 1
            
            else:
                primary = "fp"
                rp_number += 1
                if rp_number == rp_candidate_count:
                    one_done += 1
        
        elif one_done == 1:
            if primary == "fp":
                fp_number += 1
                if fp_number == fp_candidate_count:
                    one_done += 1

            else:
                rp_number += 1
                if rp_number == rp_candidate_count:
                    one_done += 1

    if one_done == 2:
        break
    
    if het_deltaG > min_dimer_value:   
        break
else:
    print("No successful primer pairs")

#If heterodimer analysis was successful, prints the final primers along with their individual statistics. Beeps when done
if het_deltaG > min_dimer_value: 
    fp_final_row = sorted_fp_candidates[sorted_fp_candidates["sequence"] == fp_hetero_anal]
    rp_final_row = sorted_rp_candidates[sorted_rp_candidates["sequence"] == rp_hetero_anal]

    print("Primers successfully generated: \nForward primer: ",fp_final_row,"\nReverse primer: ",rp_final_row,"\nHeterodimer value: ",het_deltaG)  

winsound.Beep(1000, 500)






