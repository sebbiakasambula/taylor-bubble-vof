import re
from config import REPO
s = open(REPO / 'paper' / 'paper.tex', encoding='utf-8').read()
i = s.index('\\begin{abstract}'); j = s.index('\\bibliographystyle')
for n, x in enumerate(re.split(r'(?<=[.])\s+(?=[A-Z\\])', s[i:j])):
    if re.search(r'\d\\%', x) and re.search(r'Viana|Llewellin|reference|correlation|above|below|differ|than|GCI|times', x):
        print('-', x.strip().replace('\n', ' ')[:400])
