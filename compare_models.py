import math

def read_file_to_list(filename):
    data = []
    with open(filename, 'r') as file:
        for line in file:
            # Split the line by any whitespace (handles multiple spaces/tabs)
            split_line = line.strip().split()
            
            # Try to convert each element to a number if possible
            processed_line = []
            for element in split_line:
                try:
                    # Convert to integer or float if possible
                    processed_line.append(float(element) if '.' in element else int(element))
                except ValueError:
                    # Keep as string if conversion is not possible
                    processed_line.append(element)

            data.append(processed_line)

    return data


model_names = ('constant', 'different', 'early', 'no', 'recent')

models = []
for mn in model_names:
    lines = read_file_to_list(f'./{mn}_geneflow/bestrun/darija_{mn}.bestlhoods')
    md = {}
    md['name'] = mn
    md['k'] = len(lines[0]) - 2
    md['max_obs'] = lines[1][-1]
    md['max_est'] = lines[1][-2]

    k = md['k']
    L = md['max_est']
    AIC = 2 * k - 2 * L    # no log assuming it's already a (negative) log-likelihood...

    md['AIC'] = AIC

    models.append(md)

sorted_models = sorted(models, key=lambda x: x['AIC'])

for md in sorted_models:
    print(md)