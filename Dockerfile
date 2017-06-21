FROM Python 2.7.13
MAINTAINER kjwang <wangk5@email.chop.edu>

ADD generatePep.py /

CMD [ "python", ""./generatePep.py" ]
