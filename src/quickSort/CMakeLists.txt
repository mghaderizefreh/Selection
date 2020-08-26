# quick sort library

SET(SORT_src ${SRCSORT}/quickSort.f90
)
# removed ${SRCSORT}/sortax_nested.f90

ADD_LIBRARY(${SORTLIB} STATIC ${SORT_src})
