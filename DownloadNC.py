import requests
from bs4 import BeautifulSoup
import urllib
import wget
import os
import time
from datetime import datetime, timedelta

#Choose the directory where the files should be stored
path = ''
os.chdir(path)

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        else:
            os.chdir(directory)
    except OSError:
        print("Error: Creating directory." + directory)

yesterday = datetime.now()-timedelta(days=1)   
yesterday = yesterday.strftime('%d %m %Y')
any_day = "" #Format dd mm yyyy                         
date = time.strptime(any_day, '%d %m %Y')  
#Change any_day to yesterday when suited

createFolder("./Data PWV {}".format(any_day)) 
print('Directory Data PWV {} created'.format(any_day))

os.chdir(path + '/Data PWV {}'.format(any_day)) 


#Enter the website here
url = ""
r = requests.get(url, allow_redirects=False)
soup = BeautifulSoup(r.content, 'html.parser')
TableContent = soup.select('table tr a')

Content = []
for i in range(0, len(TableContent)):
    Contentt = TableContent[i]['href'].split('/')[-1]
    Content.append(Contentt)

if len(str(date[7]))==2 or 1:
    subs = '{}'.format(date[0]) +'0' + '{}'.format(date[7])
else:
    subs = '{}'.format(date[0]) + '{}'.format(date[7]) 
FileName = [i for i in Content if subs in i]

urls = []
for i in range(0, len(FileName)):
    urlstemp = url + "/" + FileName[i] 
    urls.append(urlstemp)

for i in range(0, len(urls)):
    if not os.path.exists('./' + FileName[i]):    
        wget.download(urls[i])
print("All files downloaded!")
