FROM python:2.7.13
MAINTAINER kjwang <wangk5@email.chop.edu>

ADD generatePep.py /

RUN pip install biopython

CMD [ "python", ""./generatePep.py" ]
