from fractions import Fraction as F
from .algebra import interval_point


def cases():
    families = ((),(F(1,2),),(F(1,4),F(3,4)),(F(1,5),F(1,2),F(4,5)),
                (F(1,5),F(2,5),F(3,5),F(4,5)),
                (F(1,100),F(1,3),F(2,3),F(99,100)),
                (F(1,2),F(1,2)),(F(1,4),F(1,4),F(3,4),F(3,4)),
                (F(49,100),F(51,100)))
    rows = []
    for family in families:
        for sign in (1,-1):
            coefficients = [F(sign)]
            for root in family:
                new = [F(0)]*(len(coefficients)+1)
                for j,value in enumerate(coefficients):
                    new[j]-=root*value
                    new[j+1]+=value
                coefficients = new
            coefficients += [F(0)]*(5-len(coefficients))
            cuts = sorted({F(0),F(1),*family})
            positive = []
            for left,right in zip(cuts,cuts[1:]):
                mid = (left+right)/2
                if sum(value*mid**j for j,value in enumerate(coefficients))>0:
                    positive.append((left,right))
            point = tuple(map(F,interval_point(positive)))
            support = sum(a*b for a,b in zip(coefficients,point))
            rows.append(dict(index=len(rows),roots=list(map(str,family)),sign=sign,
                             coefficients=list(map(str,coefficients)),
                             point=list(map(str,point)),support=str(support)))
    return rows
