import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


if __name__ == '__main__':
    
    rectangles = []
    #iccad2015_input.case1
    with open('input.in', 'r') as f:
        
        for _ in range(3):
            f.readline()
        
        count = 0
        
        for line in f:
            rectangle = line[:-1].split(',')
            rectangle = [int(x) for x in rectangle]
            
            rectangles.append(rectangle)
    with open('draw_color.in', 'r') as f:
        i=0
        for line in f:
            rectangles[i].append(int(line[0]))
            rectangles[i].append(i+1)
            i+=1
    print(rectangles)  
    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.autoscale(enable=True, axis='both', tight=None)
    for i in rectangles:
        x = i[0] ; y = i[1] ; width = abs(i[0]-i[2]); height = abs(i[1]-i[3]); 
        if i[4]==0:
            rect = mpatches.Rectangle((x, y), width, height, linewidth=1, edgecolor='black', facecolor='#AAAAAA88')#gray non-color
        elif i[4]==1:
            rect = mpatches.Rectangle((x, y), width, height, linewidth=1, edgecolor='black', facecolor='#82E0AA')#green
        elif i[4]==2:
            rect = mpatches.Rectangle((x, y), width, height, linewidth=1, edgecolor='black', facecolor='#85C1E9')#blue
        ax.text(x + width / 2, y + height / 2, str(i[5]), ha='center', va='center', color='black',fontsize=4)
        patch = ax.add_patch(rect)
    #display plot
    plt.show()

            
            