import sys
############################
# Filter NCBI's id by prefix |---> Ex: To obtain top NCBI ids from Megan Blast
############################
# awk -F ' ' '{print $2}' filter_megan_blast | uniq -u > uniq_filter_megan
# Filters https://www.ncbi.nlm.nih.gov/Sequin/acc.html

def main():
	accesion_prefixes = "H,N,T,R,W,AA,AI,AW,BE,BF,BG,BI,BM,BQ,BU,CA,CB,CD,CF,CK,CN,CO,CV,CX,DN,DR,DT,DV,DW,DY,EB,EC,EE,EG,EH,EL,ES,EV,EW,EX,EY,FC,FD,FE,FF,FG,FK,FL,GD,GE,GH,GO,GR,GT,GW,HO,HS,JG,JK,JZ,U,AF,AY,DQ,EF,EU,FJ,GQ,GU,HM,HQ,JF,JN,JQ,JX,KC,KF,KJ,KM,KP,KR,KT,KU,KX,KY,MF,MG,AE,CP,CY,B,AQ,AZ,BH,BZ,CC,CE,CG,CL,CW,CZ,DU,DX,ED,EI,EJ,EK,ER,ET,FH,FI,GS,HN,HR,JJ,JM,JS,JY,KG,KO,KS,AC,DP,I,AR,DZ,EA,GC,GP,GV,GX,GY,GZ,HJ,HK,HL,KH,S,AD,AH,BC,BT"
	prefixes = accesion_prefixes.split(',')

	megan_blast = sys.argv[1]

	with open(sys.argv[2], 'w') as fout:
		with open(megan_blast, 'r') as filter_file:
			for line in filter_file:
				for prefix in prefixes:
					if prefix in line[0:len(prefix)]:
						fout.write(line)

# ---------------------

if __name__ == "__main__":
	if(len(sys.argv) != 3):
		print("ERROR | USE: python filter_megan_top.py <megan_blast_grep> <output>")
		quit()

	main()

