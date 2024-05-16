from flask import Flask, render_template
import sympy as sp

app = Flask(__name__)

@app.route('/')
def calculate():

    x, y = sp.symbols('x y')
    m11 = sp.cosh(x)
    m12 = 1
    m21 = 1
    m22 = sp.tanh(x)
    
   
    g = sp.Matrix([[m11, m12], [m21, m22]])

    gamma = christoffel_symbols(g, [x, y])


    result = "Here are the Christoffel symbols are calculated for your metric:<br><br>"
    for k in range(2):
        for i in range(2):
            for j in range(2):
                result += f'\\(\\Gamma^{{{k}}}_{{{i}{j}}} = {sp.latex(sp.simplify(gamma[k][i][j]))}\\)&nbsp;&nbsp;<br>'
    
    description = """
    The Christoffel symbols \\(\\Gamma^k_{ij}\\) are used in differential geometry to express the connection coefficients of a Levi-Civita connection.
    They are given by the formula:
    \\[
    \\Gamma^k_{ij} = \\frac{1}{2} g^{kl} \\left( \\frac{\\partial g_{lj}}{\\partial x^i} + \\frac{\\partial g_{li}}{\\partial x^j} - \\frac{\\partial g_{ij}}{\\partial x^l} \\right)
    \\]
    Here, \\(g_{ij}\\) is the metric tensor, and \\(g^{kl}\\) is its inverse. The indices \\(i\\), \\(j\\), and \\(k\\) run over the dimensions of the manifold.
    """

    return render_template('result.html', result=result, description=description)

def christoffel_symbols(g, coords):
    n = g.shape[0]
    gamma = [[[0]*n for _ in range(n)] for _ in range(n)]
    
    # Check for diagonal matrix
    if g[0, 1] == 0 and g[1, 0] == 0:
        inv_g = sp.diag(1/g[0,0], 1/g[1,1])
    else:
        inv_g = g.inv()

    for k in range(n):
        for i in range(n):
            for j in range(n):
                gamma[k][i][j] = 0.5 * sum(inv_g[k, l] * (sp.diff(g[l, j], coords[i]) + sp.diff(g[l, i], coords[j]) - sp.diff(g[i, j], coords[l])) for l in range(n))

    return gamma

if __name__ == '__main__':
    app.run(debug=True)
