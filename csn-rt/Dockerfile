# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster

WORKDIR /app

RUN pip3 install covseisnet h5py xarray dash
RUN pip3 install scipy==1.8.1

COPY csn_rt.py .
COPY run.py .

EXPOSE 8050
ENTRYPOINT [ "python", "-u", "run.py"]
CMD ["args"]