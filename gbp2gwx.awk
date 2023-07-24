#!/usr/bin/awk -f

# Converter from Generalized Bezier Patch format

NR == 1 {
    n = $1
    d = $2
    side = col = row = 0
    s_prev = n - 1
    s_next = 1
    print n
}

NR > 2 {
    if (col >= d - row) {
        if (++side >= n) {
            side = 0
            row++
        }
        s_prev = (side + n - 1) % n
        s_next = (side + 1) % n
        col = row
    }
    cp[side,col,row] = $0
    cp[s_prev,d-row,col] = $0
    cp[s_next,row,d-col] = $0
    col++
}

END {
    for (i = 0; i < n; i++) {
        for (k = 0; k < 2; k++) {
            print d
            knots = 2 * (d + 1)
            for (j = 0; j <= d; j++)
                knots = knots " 0"
            for (j = 0; j <= d; j++)
                knots = knots " 1"
            print knots
            print d + 1
            for (j = 0; j <= d; j++)
                if (k == 0)
                    print cp[i,j,k]
                else
                    print scale(cp[i,j,0], cp[i,j,1])
        }
    }
}

function scale(p, q) {
    split(p, a)
    split(q, b)
    return a[1] + (b[1] - a[1]) * d " " \
           a[2] + (b[2] - a[2]) * d " " \
           a[3] + (b[3] - a[3]) * d
}
