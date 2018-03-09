import sys
import math
from PIL import Image, ImageFont, ImageDraw

def matrix_to_png(read_mat, contig_mat, max_hits, img_dim):
	MAX_COLOR = 255
	FACTOR = float(MAX_COLOR / float(max_hits))

	img = Image.new("RGB", (img_dim, img_dim))

	draw = ImageDraw.Draw(img)
	font_large = ImageFont.load_default()

	for i in range(0, len(read_mat)-1, 1):
		curr_y = int(math.floor(i / img_dim))
		curr_x = i % img_dim
		r_val = int(math.floor(read_mat[i]*FACTOR))
		c_val = int(math.floor(contig_mat[i]*FACTOR))
		#print(str(i) + " : " + str(read_mat[i]) + "-" + str(contig_mat[i]))
		if r_val > c_val:
			img.putpixel((curr_x, curr_y), (c_val, r_val, c_val))
		else:
			img.putpixel((curr_x, curr_y), (c_val, r_val, r_val))
		
		if r_val == c_val and r_val == 0:
			img.putpixel((curr_x, curr_y), (0, 0, 0))
		#img.putpixel((curr_x, curr_y), (255-r_val, 255-c_val, 0))

	"""
	draw.rectangle(((img_dim/10 - 10, img_dim+40), (img_dim/10 + 10, img_dim+60)), fill="red")
	draw.rectangle(((3*img_dim/10 - 10, img_dim+40), (3*img_dim/10 + 10, img_dim+60)), fill="blue")
	draw.text((img_dim/10, img_dim+50), "Reads", (255,255,255), font=font_large)
	draw.text((3*img_dim/10, img_dim+50), "Contig", (255,255,255), font=font_large)
	draw.text((5*img_dim/10 + 50, img_dim+50), "Both", (255,255,255), font=font_large)
	"""

	img.save("RAC_coverage.png")

def load_seq_file(seq_fname):
	matrix_list = list()
	max_hits = 0
	with open(seq_fname, 'r') as f:
		lines = f.readlines()

		max_hits = lines[0].split(';')[-1].split('<')[-1];

		for line in lines[1:]:
			info = line.replace('\n','').split(';')
			matrix_list.append(int(info[-1]))

	return matrix_list, int(max_hits)


# -----------

############################
# Create Scaled Matrix from 2 '.seq' files (from uniseqDBCoverage) + Draw + Coverage Comparison
############################

def main():
	read_seq_fname = sys.argv[1]
	contig_seq_fname = sys.argv[2]
	dim = 1024

	r_mat, r_mhits = load_seq_file(read_seq_fname)
	c_mat, c_mhits = load_seq_file(contig_seq_fname)

	matrix_to_png(r_mat, c_mat, r_mhits, dim)

# ---------------------

if __name__ == "__main__":
	if(len(sys.argv) != 3):
		print("ERROR | USE: python drawPixels.py <1.seq> <2.seq>"); quit()

	main()