import os
import re

# set up regular token
# use https://regexper.com to visualise these if required
rx_dict = {
    'school': re.compile(r'School = (?P<school>.*)\n'),
    'grade': re.compile(r'Grade = (?P<grade>\d+)\n'),
    'name_score': re.compile(r'(?P<name_score>Name|Score)'),	
    #----------------------------------------------------------#
    'subckt_name': re.compile(r'.subckt (\S+)'),	# Capture one-or-more characters of non-whitespace after the initial match
}

def _parse_line(line):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex
    """

    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None
  
  
def parse_spice_line(line):
    """
    Do a whitespace split on the line and
    check what is going on. 
    """
    ret = []
    token = line.split()
    if len(token):
      if token[0].startswith('*'):
	print line
      
      if token[0] == ".subckt":
	subckt_name = token[1]
	NofOutConns = 0
	for outConn in token[2:]:
	  NofOutConns += 1 # count number of outConns
	ret.append(subckt_name)
	ret.append(NofOutConns)
      
    return ret
  
def parse_spice_file(filepath):
    """
    Parse text at given filepath

    Parameters
    ----------
    filepath : str
        Filepath for file_object to be parsed

    Returns
    -------
    data : pd.DataFrame
        Parsed data

    """

    data = []  # create an empty list to collect the data
    # open the file and read through it line by line
    with open(filepath, 'r') as file_object:
        line = file_object.readline()
        while line:
            # at each line check for a match with a regex
            key, match = _parse_line(line)