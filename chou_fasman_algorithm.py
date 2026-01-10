p_alpha = {
    "Glu": 1.53, "Ala" : 1.45, "Leu":1.34, "His":1.24,
    "Met": 1.20, "Gln" : 1.17, "Trp":1.14, "Val":1.14,
    "Phe": 1.12, "Lys" : 1.07, "Ile":1.00, "Asp":0.98,
    "Thr": 0.82, "Ser" : 0.79, "Arg":0.79, "Cys":0.77,
    "Asn": 0.73, "Tyr" : 0.61, "Pro":0.59, "Gly":0.53
}
p_beta = {
    "Met" : 1.67, "Val" : 1.65, "Ile":1.60, "Cys":1.30,
    "Tyr" : 1.29, "Phe" : 1.28, "Gln":1.23, "Leu":1.22,
    "Thr" : 1.20, "Trp" : 1.19, "Ala":0.97, "Arg":0.90,
    "Gly" : 0.81, "Asp" : 0.80, "Lys":0.74, "Ser":0.72,
    "His" : 0.71, "Asn" : 0.65, "Pro":0.62, "Glu":0.26
}
w_tos = { # Mapping is done from single-letter to three-letter amino acid
    "A":"Ala", "R": "Arg","N": "Asn", "D": "Asp", "C": "Cys",
    "E":"Glu", "Q": "Gln", "G": "Gly", "H":"His", "I":"Ile",
    "L": "Leu", "K":"Lys", "M": "Met", "F":"Phe", "P": "Pro",
    "S": "Ser", "T":"Thr", "W": "Trp", "Y": "Tyr","V":"Val"
}
inp_seq = ("MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYLEKPMEIARIVARCLWEESRLLQTAATAAQQGGQANHPTAAVVTEKQQMLEQHLQDVRKRVQDLEQKMKVVENLQDDFDFNYKTLKSQGDMQDLNGNNQSVTRQKMQQLEQMLTALDQMRRSIVSELAGLLSAMEYVQKTLTDEELADWKRRQQIACIGGPPNICLDRLENWITSLAESQLQTRQQIKKLEELQQKVSYKGDPIVQHRPMLEERIVELFRNLMKSAFVVERQPCMPMHPDRPLVIKTGVQFTTKVRLLVKFPELNYQLKIKVCIDKDSGDVAALRGSRKFNILGTNTKVMNMEESNNGSLSAEFKHLTLREQRCGNGGRANCDASLIVTEELHLITFETEVYHQGLKIDLETHSLPVVVISNICQMPNAWASILWYNMLTNNPKNVNFFTKPPIGTWDQVAEVLSWQFSSTTKRGLSIEQLTTLAEKLLGPGVNYSGCQITWAKFCKENMAGKGFSFWVWLDNIIDLVKKYILALWNEGYIMGFISKERERAILSTKPPGTFLLRFSESSKEGGVTFTWVEKDISGKTQIQSVEPYTKQQLNNMSFAEIIMGYKIMDATNILVSPLVYLYPDIPKEEAFGKYCRPESQEHPEADPGSAAPYLKTKFICVTPTTCSNTIDLPMSPRTLDSLMQFGNNGEGAEPSAGGQFESLTFDMELTSECATSPM")
seq = [];
for x in inp_seq:# Converting single-letter sequence to three-letter amino acid
    seq.append(w_tos[x])
n = len(seq) # length of input sequence given
helix = []
beta_strand = []
temp = ['']*n # Each index can have 'H', 'S', or '' it is a temporary array
for i in range(0,n-5): #Helix prediction code
    hc= 0
    for j in range(i, i+6):
        if (p_alpha[seq[j]]) >= 1.00:
            hc += 1
    if hc >= 4: #extend if more than or equal to 4(avg prospenity), have prospenty greater than or equal to 1. of 4 out 6 residues
        start = i
        end = i + 6 
        while start > 0 and start+2<n: #left extension
            remaning_part = seq[start-1:start+3]
            avga = sum(p_alpha[x] for x in remaning_part) / 4
            if avga >= 1.0:
                start -= 1
            else:
                break
        while end+1<=n and end-3>=0: # right extension
            remaning_part = seq[end-3:end+1]
            avga = sum(p_alpha[x] for x in remaning_part) / 4
            if avga >= 1.0:
                end += 1
            else:
                break
        helix.append((start, end-1))
        for k in range(start, end):
            if 'H' not in temp[k]:
                temp[k] +='H' # stored regions which are predicted to be helix

for i in range(0, n-4): # beta strand prediction
    b_c = 0
    for j in range(i, i+5):
        if p_beta[(seq[j])] >= 1.00:
            b_c += 1

    if b_c>= 3: #extend if more than or equal to 3(avg prospenity), have prospenty greater than or equal to 1. of 3 out 5 residues
        start = i
        end = i + 5
        while start > 0 and start+2<n: #left extension, 3 from old and one new amino acid towards left
            remaning_part = seq[start-1:start+3]
            avgb = sum(p_beta[x] for x in remaning_part) / 4
            if avgb >= 1.00:
                start -= 1
            else:
                break
        while end+1 <= n and end-3>=0: # right extension, 3 from old and one new amino acid towards right
            remaning_part = seq[end-3:end+1]
            avgb = sum(p_beta[x] for x in remaning_part) / 4
            if avgb >= 1.00:
                end += 1
            else:
                break
        beta_strand.append((start, end-1))
        for k in range(start, end):
            if 'S' not in temp[k]:
                temp[k] += 'S' # stored predicted to be beta strands
#all helix and strand regions BEFORE conflict resolution (including conflicts)
ph_r = []
pb_r = []
start= -1
for i in range(n): # helix including conflict
    if 'H' in temp[i]:
        if start == -1:
            start= i
    else:
        if start!= -1:
            ph_r.append((start, i-1))
            start= -1
if start!= -1:
    ph_r.append((start, n-1))

start =-1
for i in range(n): # beta strand including conflict
    if 'S' in temp[i]:
        if start== -1:
            start= i
    else:
        if start!= -1:
            pb_r.append((start, i-1))
            start = -1
if start!= -1:
    pb_r.append((start, n-1))

conflict_pos =[] # conflict resolving and detection is done here
cp_seg =[] # it store segments like from (0,22) are conflicted region
resolve_info= [] # it stores resolved conflicts
for i in range(len(temp)):
    val= temp[i]
    if 'H' in val and 'S' in val: # position in temp which are conflicted denoted by HS is identified here
        conflict_pos.append(i)
if conflict_pos: # grouping of consecutive conflicting position is done here then stored in cpseg
    sstart =conflict_pos[0]
    prev =sstart
    for i in range(1, len(conflict_pos)):
        pos = conflict_pos[i]
        if pos == prev + 1:  
            prev = pos      
        else:
            cp_seg.append((sstart, prev))  
            sstart = pos 
            prev = pos
    cp_seg.append((sstart, prev))
temp2 =['-'] * n   # conflict resolving, one more temporary array is created
for i in range(n):
    if temp[i]== 'H':
        temp2[i] = 'H'
    elif temp[i]== 'S':
        temp2[i] ='S'
  
for (a, b) in cp_seg:
    block= seq[a:b+1]
    s_pa =sum(p_alpha[x] for x in block) # here it check conflicted region, by comparing sum of prospenity and assigning higher propority based on higher value
    s_pb =sum(p_beta[x] for x in block)
    if s_pa>=s_pb:
        c = 'H'
    else:
        c= 'S'
    for k in range(a, b+1):
        temp2[k] = c
    resolve_info.append((a, b, s_pa, s_pb, c))

for i in range(n):
    if temp2[i] not in ('H', 'S'): # remaning position which are nothing(turn etc) denoted by -
        temp2[i]= '-'
f_helix =[] # final helix region is stored
start = -1
for i in range(n):
    if temp2[i] == 'H':
        if start== -1:
            start= i
    else:
        if start != -1:
            f_helix.append((start,i - 1))
            start = -1
if start!= -1:
    f_helix.append((start, n - 1))

f_beta = [] # final beta strand region is stored
start= -1
for i in range(n):
    if temp2[i] == 'S':
        if start == -1:
            start = i
    else:
        if start != -1:
            f_beta.append((start,i - 1))
            start = -1
if start!= -1:
    f_beta.append((start,n - 1))
# printing of results and statistics ## I have taken 0 based array indexing
print("a) Helix regions (including conflicted regions):")
if ph_r:
    for (a, b) in ph_r:
        print(f"Helix: positions ({a}-{b}) -> {''.join(seq[a:b+1])}") # these are position that are predicted to be helix. it contains conflicted region as well
else:
    print("No helix regions found.")
print()

print("b) Beta strand regions (including conflicted regions):") 
if pb_r:
    for (a, b) in pb_r:
        print(f"Beta strand: positions ({a}-{b}) -> {''.join(seq[a:b+1])}") # these are position that are predicted to be beta strand. it contains conflicted region as well
else:
    print("No beta strand regions found.")
print()

print("Helical regions (final after conflict resolution):")
if f_helix:
    for(g,z) in f_helix:
        print(f"Helix: positions ({g}-{z})-> {''.join(seq[g:z+1])}") #these are final helix position
else:
    print("No helix region is found.")
print()

print("Beta strand regions (final after conflict resolution):")
if f_beta:
    for(g, z) in f_beta:
        print(f"Beta strand: positions ({g}-{z}) ->{''.join(seq[g:z+1])}") #these are final beta strand position
else:
    print("No beta strand region is found ")
print()
print("c)Conflicting region:")
if cp_seg:
    for (a, b, s_pa, s_pb, c) in resolve_info:
        print(f"Conflict positions({a}-{b}) -> seq = {''.join(seq[a:b+1])}") # conflicted regions
        print(f"sum P_alpha = {s_pa: }, sum P_beta = {s_pb: }  => assigned to {c}") # conflict resolution
else:
    print("No conflicts detected")
print()
print("\nFinal full-sequence allignment:\n")
print("Sequence : ", ''.join(seq))
print("Structure: ", ''.join(temp2))