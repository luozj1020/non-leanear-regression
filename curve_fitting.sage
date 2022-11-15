# for this example, we first generate some data with built-in variance:
x = [-1, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
y = [-0.8669, -0.2997,  0.3430, 1.0072, 1.6416, 2.2022, 2.6558, 2.9823, 3.1755, 3.2416, 3.1974]
data = [(x[i], y[i]) for i in range(len(x))]

# design a model with adjustable parameters a,b,c that describes the data
var('a, b, c, d, e, x')
model(x) = a + b * sin(c * x) + d * cos(e * x)

# calculate the values of the parameters that best fit the model to the data
sol = find_fit(data,model)
show(sol)

# define f(x), the model with the parameters set to the fitted values
f(x) = model(a=sol[0].rhs(),b=sol[1].rhs(),c=sol[2].rhs(),d=sol[3].rhs(),e=sol[4].rhs())

# create an empty plot object
a = plot([])
# add a plot of the model, with respect to x from -pi/2 to pi/2
a += plot(f(x),x,[-pi/2,pi/2])

# add a plot of the data, in red
a += list_plot(data,color='red')
show(a)

e = 0
for i in range(len(x)):
    e += (f(x[i])-y[i])**2
print(e)