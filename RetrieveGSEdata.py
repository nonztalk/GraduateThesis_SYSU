import urllib2

# example: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12435&targ=self&form=text&view=full
def get_GSE_address(GSE_No):
    url_base = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
    access = GSE_No # GSE number to require
    scope = "targ=self"
    form = "form=text"
    text = "view=full"

    url = url_base + access + '&' + scope + '&' + form + '&' + text

    f = urllib2.urlopen(url)
    data = f.read()

    with open("GSE_download.txt", "wb") as dataFile:
        for line in data:
            if "!Series_supplementary_file" in line:
                address = line.lstrip("!Series_supplementary_file = ")
                dataFile.write(string.replace("ftp://ftp.ncbi.nlm.nih.gov", "anonftp@ftp-private.ncbi.nlm.nih.gov"))
                dataFile.write("\n")
