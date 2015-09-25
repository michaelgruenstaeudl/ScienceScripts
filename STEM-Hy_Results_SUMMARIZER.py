import re
import sys
p="/home/michael_science/git/github.com_ScienceScripts/"
sys.path.insert(0, p)
print "\nLibrary path:", p, "\n"
import GeneralFileOperations as GFO
import GeneralStringOperations as GSO

inFn = sys.argv[1]
inData = GFO.loadR(inFn)

handle = inData.split(">")
# Remove all empty strings
handle = filter(None, handle)
# Strip all line breaks from list elements
handle = [elem.rstrip("\n") for elem in handle]

out_list = []
for elem in handle:
    # indcs tell where numbers begin
    indcs = [m.start()+5 for m in re.finditer('AIC: ', elem)]
    # extract numbers
    values = [elem[i:i+10] for i in indcs]
    # Convert a list of strings into a list of integers
    values = map(float, values)
    # Ensure that all list values negative
    for i,v in enumerate(values):
        if v > 0:
            values[i] = -abs(v)

    out_handle = elem.split("\n")[0]

    if len(values) == 0:
        print("Error: " + out_handle)
    else:
        keyval = values.index(min(values))
        tmp_data = out_handle.split("_")
        # DEBUG  print "\n", tmp_data, "\n"
        if keyval == 2:
            out_list.append(",".join([tmp_data[1]] + [tmp_data[0]] + ["1"]))
        if keyval < 2:
            out_list.append(",".join([tmp_data[1]] + [tmp_data[0]] + ["0"]))

outData = '\n'.join(out_list)

tmpFn = GSO.rmext(inFn)+".hybridplot"
GFO.save(tmpFn, outData)

