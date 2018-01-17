FROM alpine
MAINTAINER <http://o2r.info>

ENV HOST o2r.info
ENV COUNT 30
ENV SECONDS 60

RUN chmod u+s /bin/ping

# ping 30 times, but at most 60 seconds
CMD ping -c $COUNT -w $SECONDS $HOST

