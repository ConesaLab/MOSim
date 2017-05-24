
is.declared <- function(object, key = NULL) {

    if (! is.null(key)) {

        if (! isS4(object)) {
            if (! exists(key, where = object))
                return(FALSE)

            object <- get(key, object)
        } else {
            object <- slot(object, key)
        }
    }

    # Check length, is.null, etc
    check <- ! is.null(object)

    # Check empty matrix
    if (check) {
        if (! is.null(dim(object))) {
            check <- ! all(is.na(object))
        } else {
            check <- as.logical(length(object))
        }
    }

    return(check)
}
