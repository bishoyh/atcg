# ftpkeggdre.py, get zebrafish xml kegg pathways

STOREDIR = 'drexml/'

from ftplib import FTP
#import os



f = FTP('ftp.genome.jp')
print "Welcome:", f.getwelcome()
f.login()

f.cwd('/pub/kegg/xml/organisms/dre')
print "CWD: ", f.pwd()

items = f.nlst()

print "%d items:" %len(items)
#print items[:7]

for ii in items:
    if '.xml' in ii:
        print ii
        localwrite = open(STOREDIR+ii, 'wb')    #use binary mode
        f.retrbinary('RETR ' + ii, localwrite.write)
        localwrite.close()



f.quit()
