


def read_sample_file(sample, circles, infile):
    I = open(infile)
    H = I.readline()
    for line in I:
	L = line.replace('\n', '').split('\t')
	if not L[0] in circles:
	    circles[L[0]] = {'transcript' : L[1]}
	circles[L[0]][sample] = [L[5], L[6], L[8], L[9]]    
    I.close()
    return(circles)

def iterate_over_sample(sample, folder):
    circles = {}
    samples = []
    for lola in sample:
	circles = read_sample_file(lola.split('.')[0], circles, '%s/%s' %(folder, lola))
	samples += [lola.split('.')[0]]
    return(circles, samples)

def write_tables(outfolder, samples,circles):
    O_s = open('%s/single_breakpoints.txt' %(outfolder), 'w')
    O_d = open('%s/double_breakpoints.txt' %(outfolder), 'w')
    O_l = open('%s/circle_annotated_length.txt' %(outfolder), 'w')
    O_c = open('%s/circle_bases_covered.txt' %(outfolder), 'w')
    sorted_samples = sorted(samples)
    O_s.write('circle_id\ttranscript\t%s\n' %('\t'.join(sorted_samples)))
    O_d.write('circle_id\ttranscript\t%s\n' %('\t'.join(sorted_samples)))
    O_l.write('circle_id\ttranscript\t%s\n' %('\t'.join(sorted_samples)))
    O_c.write('circle_id\ttranscript\t%s\n' %('\t'.join(sorted_samples)))
    sorted_circles = sorted(circles.keys())
    for lola in sorted_circles:
	O_s.write('%s\t%s' %(lola, circles[lola]['transcript']))
	O_d.write('%s\t%s' %(lola, circles[lola]['transcript']))
	O_l.write('%s\t%s' %(lola, circles[lola]['transcript']))
	O_c.write('%s\t%s' %(lola, circles[lola]['transcript']))
	for s in sorted_samples:
	    if s in circles[lola]:
		O_s.write('\t%s' %circles[lola][s][0])
		O_d.write('\t%s' %circles[lola][s][1])
		O_l.write('\t%s' %circles[lola][s][2])
		O_c.write('\t%s' %circles[lola][s][3])
	    if not s in circles[lola]:
		O_s.write('\t0')
		O_d.write('\t0')
		O_l.write('\t0')
		O_c.write('\t0')
	O_s.write('\n')
	O_d.write('\n')
	O_l.write('\n')
	O_c.write('\n')
    O_s.close()
    O_d.close()
    O_l.close()
    O_c.close()
    return


frontal_cortex = ('frontal_cortex_A.mate_status.added.txt', 'frontal_cortex_B.mate_status.added.txt', 'frontal_cortex_C.mate_status.added.txt', 'frontal_cortex_D.mate_status.added.txt', 'frontal_cortex_E.mate_status.added.txt', 'frontal_cortex_F.mate_status.added.txt')
mouse = ('old_cerebellum.mate_status.added.txt', 'old_hippocampus.mate_status.added.txt', 'old_liver.mate_status.added.txt', 'young_cerebellum.mate_status.added.txt', 'young_hippocampus.mate_status.added.txt', 'young_liver.mate_status.added.txt')
mouse_folder = '/home/fmetge/Documents/work/circRNA/FUCHS/ECCB/mouse/FUCHS'
frontal_cortex_folder = '/home/fmetge/Documents/work/circRNA/FUCHS/ECCB/mouse/frontal_cortex'


FC, samples_FC = iterate_over_sample(frontal_cortex, frontal_cortex_folder)
M, samples_M = iterate_over_sample(mouse, mouse_folder)


write_tables(frontal_cortex_folder, samples_FC, FC)
write_tables(mouse_folder, samples_M, M)
