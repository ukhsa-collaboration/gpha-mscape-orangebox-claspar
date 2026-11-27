FROM python:3.14

ADD . /app
WORKDIR /app
RUN pip install .

CMD ["claspar"]
