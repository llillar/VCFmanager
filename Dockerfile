FROM python:3.7-slim
CMD ["/bin/bash"]

ENV PYTHONUNBUFFERED 1

RUN mkdir /code
WORKDIR /code

COPY requirements.txt /code
RUN python3.7 -m pip install --upgrade pip \
&& python3.7 -m pip install -r requirements.txt
COPY /code/ /code/