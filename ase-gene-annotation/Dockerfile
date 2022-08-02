FROM python

LABEL org.opencontainers.image.source https://github.com/icgc-argo-workflows/allele-specific-expression

ENV PATH="/tools:${PATH}"

COPY *.py /tools/

ENTRYPOINT ["/usr/bin/env"]

CMD ["/bin/bash"]
