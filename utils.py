def get_params_dict(param_data):
    params = {}
    for entry in param_data:
        if entry[1] == 'False':
            params[entry[0]] = False
        elif entry[1] == 'True':
            params[entry[0]] = True
        else:
            try:
                params[entry[0]] = int(entry[1])
            except:
                params[entry[0]] = entry[1]    
    
    return params

def calculate_padding(inputLength, filterSize):
    padding = inputLength - (inputLength - filterSize + 1)
    return int(padding/2) #appended to both sides the half of what is needed