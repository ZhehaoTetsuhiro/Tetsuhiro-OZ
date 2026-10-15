package optics

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Expression language for user-defined media (自定义介质表达式)
//
// A custom medium (see medium_custom.go) describes its refractive index
// n(x, y, z, λ) — and optionally an absorption/gain coefficient k(x, y, z, λ) —
// as an arithmetic expression. This file compiles such an expression onto the
// medium's own evaluation environment, which carries a z coordinate that the
// thin-element environment of expr.go does not have.
//
// The language itself is shared with the scripted-element language of expr.go:
// the lexer, the recursive-descent parser (parseExprNode), the operator set,
// the function table (exprFuncs) and ExprError with its source columns all come
// from there. Only the *compiler* differs — it turns the AST into a closure over
// mediumEnv instead of exprEnv — so that a medium can read the depth z and the
// per-layer depth zl. The grammar is identical:
//
//	expr    := or
//	or      := and ('||' and)*
//	and     := cmp ('&&' cmp)*
//	cmp     := add (('>'|'<'|'>='|'<='|'=='|'!=') add)*
//	add     := mul (('+'|'-') mul)*
//	mul     := unary (('*'|'/'|'%') unary)*
//	unary   := ('-'|'+'|'!') unary | power
//	power   := primary ('^' unary)?          // right-associative
//	primary := number | ident | ident '(' args ')' | '(' expr ')'
//
// Built-in names: x, y (transverse coordinates, m, in the grid frame), z (depth
// from the medium's entrance face, m), zl (depth from the current layer's
// entrance face, m), r (=hypot(x,y)), th (=atan2(y,x), rad), wl (wavelength, m),
// k (=2*pi/wl), pi, e. Any other identifier must be bound by the layer's
// "params" map as a compile-time constant; an unknown name is a load-time
// error, never a silent zero.
//
// The language has no loops, no arrays, no assignment and no I/O: an expression
// cannot hang and cannot read anything but its own inputs, so evaluating it
// twice at the same point gives the same value. That is what makes it safe to
// run user-supplied text inside the propagation loop: the expression is parsed
// and compiled *once* (at element-construction time) and thereafter evaluated as
// a tree of plain float64 closures, never as a string.
// ---------------------------------------------------------------------------

// mediumEnv is the per-point evaluation environment of a medium expression. r
// and th are derived from x,y only when the expression asks for them.
//
// frames is the argument arena used by the function-call nodes: the expression
// compiler gives every call node its own three-float slot index (call nodes are
// never nested in the same slot, so a nested call cannot clobber its parent's
// arguments), and the arena itself belongs to the caller, i.e. to one
// propagation. That keeps the per-pixel path allocation-free without making a
// compiled medium unsafe to use from two goroutines at once.
type mediumEnv struct {
	x, y, z, zl, r, th, wl, k float64
	frames                    []float64
}

// mediumFrameWidth is the number of argument slots reserved per call node (3 =
// the largest arity in exprFuncs).
const mediumFrameWidth = 3

type mediumFunc func(*mediumEnv) float64

// MediumExpr is a parsed and compiled medium expression.
type MediumExpr struct {
	src    string
	root   mediumFunc
	needsR bool
	needsT bool
	// spatial is true when the expression reads x, y, z, zl, r or th, i.e. when
	// its value can vary from point to point inside a layer. A non-spatial
	// expression (constants, params and the wavelength only) describes a
	// homogeneous slice, which the split-step can propagate exactly.
	spatial bool
	names   []string
	// nFrames is the number of function-call nodes in the tree, i.e. the
	// number of mediumFrameWidth-sized slots an environment must provide.
	nFrames int
}

// FrameSlots returns the number of float64 slots an environment needs to
// evaluate the expression (MediumFrameSlots * nFrames).
func (m *MediumExpr) FrameSlots() int { return m.nFrames * mediumFrameWidth }

// newMediumEnv builds an evaluation environment with an arena large enough for
// the expression. Callers that evaluate many expressions can size one arena
// themselves and reuse it (see PropagateCustomMedium).
func (m *MediumExpr) newMediumEnv(wl float64) *mediumEnv {
	e := &mediumEnv{wl: wl, k: 2 * math.Pi / wl}
	if m.nFrames > 0 {
		e.frames = make([]float64, m.FrameSlots())
	}
	return e
}

// Source returns the expression text it was compiled from.
func (m *MediumExpr) Source() string { return m.src }

// Names lists the variables and parameters the expression references, sorted.
func (m *MediumExpr) Names() []string { return append([]string(nil), m.names...) }

// IsSpatial reports whether the value can vary in space (it reads at least one
// of x, y, z, zl, r, th).
func (m *MediumExpr) IsSpatial() bool { return m.spatial }

// NeedsRadius reports whether the expression reads r.
func (m *MediumExpr) NeedsRadius() bool { return m.needsR }

// NeedsTheta reports whether the expression reads th.
func (m *MediumExpr) NeedsTheta() bool { return m.needsT }

// Eval evaluates the expression at one point. z is the depth from the medium's
// entrance face and zl the depth from the current layer's entrance face (both
// in metres); wl is the wavelength in metres.
func (m *MediumExpr) Eval(x, y, z, zl, wl float64) float64 {
	e := m.newMediumEnv(wl)
	e.x, e.y, e.z, e.zl = x, y, z, zl
	return m.evalEnv(e)
}

// evalEnv evaluates against a caller-managed environment, so the hot loop
// reuses one struct and one set of coordinates instead of allocating per point.
// A caller that does not size the argument arena itself gets one here (the hot
// loop pre-sizes it once for the widest layer expression, so the check below
// passes immediately).
func (m *MediumExpr) evalEnv(e *mediumEnv) float64 {
	if m.nFrames > 0 && len(e.frames) < m.FrameSlots() {
		e.frames = make([]float64, m.FrameSlots())
	}
	if m.needsR || m.needsT {
		e.r = math.Hypot(e.x, e.y)
		if m.needsT {
			e.th = math.Atan2(e.y, e.x)
		}
	}
	return m.root(e)
}

// mediumBuiltins are the variables every medium expression may use.
var mediumBuiltins = map[string]mediumFunc{
	"x":  func(e *mediumEnv) float64 { return e.x },
	"y":  func(e *mediumEnv) float64 { return e.y },
	"z":  func(e *mediumEnv) float64 { return e.z },
	"zl": func(e *mediumEnv) float64 { return e.zl },
	"r":  func(e *mediumEnv) float64 { return e.r },
	"th": func(e *mediumEnv) float64 { return e.th },
	"wl": func(e *mediumEnv) float64 { return e.wl },
	"k":  func(e *mediumEnv) float64 { return e.k },
	"pi": func(*mediumEnv) float64 { return math.Pi },
	"e":  func(*mediumEnv) float64 { return math.E },
}

// mediumBuiltinNames lists the built-in variable names (documentation, catalog).
func mediumBuiltinNames() []string {
	return []string{"x", "y", "z", "zl", "r", "th", "wl", "k", "pi", "e"}
}

// MediumExprBuiltins returns the built-in variable names of the medium
// expression language, sorted, for documentation and error messages.
func MediumExprBuiltins() []string { return mediumBuiltinNames() }

// MediumExprFunctions returns the callable function names shared by the medium
// and scripted-element expression languages (expr.go's function table).
func MediumExprFunctions() []string {
	out := make([]string, 0, len(exprFuncs))
	for name := range exprFuncs {
		out = append(out, name)
	}
	sortStrings(out)
	return out
}

// CompileMediumExpr parses src and compiles it, binding every non-built-in
// identifier to its constant value in params. The result is evaluated as
// closures, so a compiled expression costs no parsing per point.
func CompileMediumExpr(src string, params map[string]float64) (*MediumExpr, error) {
	node, err := parseExprNode(src)
	if err != nil {
		return nil, err
	}
	vars := make(map[string]mediumFunc, len(params))
	for k, v := range params {
		val := v
		vars[k] = func(*mediumEnv) float64 { return val }
	}
	var names []string
	spatial := false
	nFrames := 0
	root, err := compileMediumNode(node, vars, &names, &spatial, &nFrames)
	if err != nil {
		return nil, err
	}
	m := &MediumExpr{src: src, root: root, spatial: spatial, nFrames: nFrames}
	for _, n := range names {
		switch n {
		case "r":
			m.needsR = true
		case "th":
			m.needsT = true
		}
	}
	sortStrings(names)
	m.names = uniqueStrings(names)
	return m, nil
}

// compileMediumNode turns an AST into a closure over mediumEnv, checking every
// name as it goes (an unknown name is an error at compile time, not a zero).
func compileMediumNode(n *exprNode, vars map[string]mediumFunc, names *[]string, spatial *bool, nFrames *int) (mediumFunc, error) {
	switch n.kind {
	case exprNum:
		v := n.num
		return func(*mediumEnv) float64 { return v }, nil
	case exprVar:
		if f, ok := vars[n.name]; ok {
			*names = append(*names, n.name)
			return f, nil
		}
		if f, ok := mediumBuiltins[n.name]; ok {
			*names = append(*names, n.name)
			switch n.name {
			case "x", "y", "z", "zl", "r", "th":
				*spatial = true
			}
			return f, nil
		}
		return nil, &ExprError{Pos: n.pos, Msg: fmt.Sprintf(
			"unknown variable %q (built-ins: %s; any other name must be given in \"params\")",
			n.name, joinNames(mediumBuiltinNames())), Name: n.name}
	case exprUn:
		a, err := compileMediumNode(n.args[0], vars, names, spatial, nFrames)
		if err != nil {
			return nil, err
		}
		switch n.op {
		case "-":
			return func(e *mediumEnv) float64 { return -a(e) }, nil
		case "!":
			return func(e *mediumEnv) float64 {
				if a(e) == 0 {
					return 1
				}
				return 0
			}, nil
		}
		return nil, &ExprError{Pos: n.pos, Msg: "unsupported unary operator " + n.op}
	case exprBin:
		a, err := compileMediumNode(n.args[0], vars, names, spatial, nFrames)
		if err != nil {
			return nil, err
		}
		b, err := compileMediumNode(n.args[1], vars, names, spatial, nFrames)
		if err != nil {
			return nil, err
		}
		switch n.op {
		case "+":
			return func(e *mediumEnv) float64 { return a(e) + b(e) }, nil
		case "-":
			return func(e *mediumEnv) float64 { return a(e) - b(e) }, nil
		case "*":
			return func(e *mediumEnv) float64 { return a(e) * b(e) }, nil
		case "/":
			return func(e *mediumEnv) float64 { return a(e) / b(e) }, nil
		case "%":
			return func(e *mediumEnv) float64 { return math.Mod(a(e), b(e)) }, nil
		case "^":
			return func(e *mediumEnv) float64 { return math.Pow(a(e), b(e)) }, nil
		case ">":
			return func(e *mediumEnv) float64 { return boolf(a(e) > b(e)) }, nil
		case "<":
			return func(e *mediumEnv) float64 { return boolf(a(e) < b(e)) }, nil
		case ">=":
			return func(e *mediumEnv) float64 { return boolf(a(e) >= b(e)) }, nil
		case "<=":
			return func(e *mediumEnv) float64 { return boolf(a(e) <= b(e)) }, nil
		case "==":
			return func(e *mediumEnv) float64 { return boolf(a(e) == b(e)) }, nil
		case "!=":
			return func(e *mediumEnv) float64 { return boolf(a(e) != b(e)) }, nil
		case "&&":
			return func(e *mediumEnv) float64 { return boolf(a(e) != 0 && b(e) != 0) }, nil
		case "||":
			return func(e *mediumEnv) float64 { return boolf(a(e) != 0 || b(e) != 0) }, nil
		}
		return nil, &ExprError{Pos: n.pos, Msg: "unsupported operator " + n.op}
	case exprCall:
		spec, ok := exprFuncs[n.name]
		if !ok {
			return nil, &ExprError{Pos: n.pos, Msg: fmt.Sprintf(
				"unknown function %q (see MediumExprFunctions for the available set)", n.name)}
		}
		if len(n.args) != spec.arity {
			return nil, &ExprError{Pos: n.pos, Msg: fmt.Sprintf("%s takes %d argument(s), got %d", n.name, spec.arity, len(n.args))}
		}
		if n.name == "if" {
			// Lazy, exactly as in expr.go: only the taken branch is evaluated,
			// so if() can guard a division or a sqrt away from its singular
			// point (e.g. if(r>0, n0+a*sin(b*r)/r, n0)).
			c, err := compileMediumNode(n.args[0], vars, names, spatial, nFrames)
			if err != nil {
				return nil, err
			}
			t, err := compileMediumNode(n.args[1], vars, names, spatial, nFrames)
			if err != nil {
				return nil, err
			}
			f, err := compileMediumNode(n.args[2], vars, names, spatial, nFrames)
			if err != nil {
				return nil, err
			}
			return func(e *mediumEnv) float64 {
				if c(e) != 0 {
					return t(e)
				}
				return f(e)
			}, nil
		}
		args := make([]mediumFunc, len(n.args))
		for i, a := range n.args {
			af, err := compileMediumNode(a, vars, names, spatial, nFrames)
			if err != nil {
				return nil, err
			}
			args[i] = af
		}
		call := spec.call
		// One arena slot per call node: nested calls use different slots, so no
		// argument buffer is ever shared between a call and its arguments. The
		// slice literal the naive version needs here would escape to the heap on
		// every pixel (measured: 3 allocations per pixel for
		// "1.5+sin(1e4*x)*cos(1e4*y)"); the arena keeps the loop allocation-free.
		frame := *nFrames * mediumFrameWidth
		*nFrames++
		n := len(args)
		return func(e *mediumEnv) float64 {
			buf := e.frames[frame : frame+n : frame+n]
			for i := range args {
				buf[i] = args[i](e)
			}
			return call(buf)
		}, nil
	}
	return nil, &ExprError{Pos: n.pos, Msg: "internal: bad expression node"}
}

// joinNames renders a name list as "a, b, c" for error messages.
func joinNames(names []string) string {
	out := ""
	for i, n := range names {
		if i > 0 {
			out += ", "
		}
		out += n
	}
	return out
}
