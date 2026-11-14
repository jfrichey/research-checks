"""Independent finite-configuration checks of Euler and boundary identities."""
from itertools import product

D4=((1,0),(-1,0),(0,1),(0,-1))
D8=tuple((a,b) for a,b in product((-1,0,1),repeat=2) if a or b)


def components(points,directions=D4):
    remaining=set(points)
    answer=[]
    while remaining:
        stack=[remaining.pop()]
        component=set(stack)
        while stack:
            x,y=stack.pop()
            for dx,dy in directions:
                v=(x+dx,y+dy)
                if v in remaining:
                    remaining.remove(v)
                    component.add(v)
                    stack.append(v)
        answer.append(component)
    return answer


def count(points):
    return len(components(points))


def main():
    euler_checks=0
    for m in range(1,5):
        for n in range(1,5):
            vertices=[(x,y) for x in range(n) for y in range(m)]
            boundary={v for v in vertices if v[0] in (0,n-1) or v[1] in (0,m-1)}
            for mask in range(1<<len(vertices)):
                black={v for i,v in enumerate(vertices) if mask>>i&1}
                white=set(vertices)-black
                edges=sum((x+dx,y+dy) in black for x,y in black for dx,dy in ((1,0),(0,1)))
                faces=sum(all((x+dx,y+dy) in black for dx,dy in ((0,0),(1,0),(0,1),(1,1)))
                          for x in range(n-1) for y in range(m-1))
                interior_white=sum(not (c & boundary) for c in components(white,D8))
                assert count(black)-interior_white == len(black)-edges+faces
                euler_checks+=1
    sandwich_checks=0
    m,n=3,4
    vertices=[(x,y) for x in range(n) for y in range(m)]
    for mask in range(1<<len(vertices)):
        black={v for i,v in enumerate(vertices) if mask>>i&1}
        actual=count(black)-count({v for v in black if v[0]<n-1})
        for k in (2,3,4):
            observed={v for v in black if v[0]>=n-k}
            old={v for v in observed if v[0]<n-1}
            wire={(n-k-1,y) for y in range(m)}
            lower=count(observed)-count(old)
            upper=count(observed|wire)-count(old|wire)
            assert lower<=actual<=upper
            sandwich_checks+=1
    print(f'Checked Euler identity on {euler_checks} finite configurations.')
    print(f'Checked {sandwich_checks} free/wired increment sandwiches.')


if __name__=='__main__':
    main()
