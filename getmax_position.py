import math
def getmax_po(fit):
	max_fit=max(fit)
	for i in range(len(fit)):
		if fit[i] == max_fit:
			return i
			break
if __name__=="__main__":
	c=[1.0,2.0,3.0,4.0,1.0,2.0]
	t=getmax_po(c)
	print(t)