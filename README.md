# Leetcode DSA in Java

## EASY

### 1. Roman to Integer - https://leetcode.com/problems/roman-to-integer
```java
public int romanToInt(String s) {
    int nums[]=new int[s.length()];
    for(int i=0;i<s.length();i++){
        switch (s.charAt(i)){
            case 'M':
                nums[i]=1000;
                break;
            case 'D':
                nums[i]=500;
                break;
            case 'C':
                nums[i]=100;
                break;
            case 'L':
                nums[i]=50;
                break;
            case 'X' :
                nums[i]=10;
                break;
            case 'V':
                nums[i]=5;
                break;
            case 'I':
                nums[i]=1;
                break;
        }
    }
    int sum=0;
    for(int i=0;i<nums.length-1;i++){
        if(nums[i]<nums[i+1])
            sum-=nums[i];
        else
            sum+=nums[i];
    }
    return sum+nums[nums.length-1];
}
```

### 2. Delete Node in Linked List - https://leetcode.com/problems/delete-node-in-a-linked-list

Linked list:
```java
public class Solution {
    public void deleteNode(ListNode node) {
		if (node.next != null) {
			node.val = node.next.val;
			node.next = node.next.next;
		}
	}
}
```

### 3. Reverse Linked List - https://leetcode.com/problems/reverse-linked-list/description/

Linked list:
```java
/**
 * Definition for singly-linked list.
 * public class ListNode {
 *     int val;
 *     ListNode next;
 *     ListNode() {}
 *     ListNode(int val) { this.val = val; }
 *     ListNode(int val, ListNode next) { this.val = val; this.next = next; }
 * }
 */
class Solution {
    public ListNode reverseList(ListNode head) {
        if(head == null) return head;

        ListNode previousNode = null;
        ListNode currentNode = head;

        while(currentNode.next != null) {
            ListNode node = currentNode.next; // save next node of next node
            currentNode.next = previousNode; // reverse pointer
            previousNode = currentNode;
            currentNode = node;
        }

        currentNode.next = previousNode;

        return currentNode;
    }
}
```

### 4. Binary Search - https://leetcode.com/problems/binary-search

Binary search:
```java
class Solution {
    public int search(int[] nums, int target) {
            return searchFunc(nums, 0, nums.length-1, target);
    }

    public int searchFunc(int[] nums, int startIndex, int endIndex, int target) {
        if(startIndex > endIndex) {
            return -1;
        }
        int mid = startIndex + (endIndex - startIndex) / 2;
        if (target == nums[mid]) {
            return mid;
        } else if (target < nums[mid]) {
            return searchFunc(nums, startIndex, mid - 1, target);
        } else {
            return searchFunc(nums, mid + 1, endIndex, target);
        }
    }
}
```

### 5. Arrange Coins in Stairs - https://leetcode.com/problems/arranging-coins/description/

Recursion:
```java
class Solution {
    public int arrangeCoins(int n) {
        return generateLevels(n, 1);
    }

    private int generateLevels(int numCoinsLeft, int currentLevel) {
        int coinsNeeded = currentLevel;

        if(numCoinsLeft < coinsNeeded) {
            // We cannot fulfill this level
            return currentLevel - 1;
        }

        return generateLevels(numCoinsLeft - coinsNeeded, currentLevel + 1);
    }
}
```

### 6. Longest Substring Without Repeating Characters - https://leetcode.com/problems/longest-substring-without-repeating-characters

Slow and fast pointer solution:
```java
public int lengthOfLongestSubstring(String s) {
    int i = 0, j = 0, max = 0;
    Set<Character> set = new HashSet<>();
    
    while (j < s.length()) {
        if (!set.contains(s.charAt(j))) {
            set.add(s.charAt(j++));
            max = Math.max(max, set.size());
        } else {
            set.remove(s.charAt(i++));
        }
    }
    
    return max;
}
```

Sliding window solution:
```java
public int lengthOfLongestSubstring(String s) {
    int left = 0; // start of window
    int right = 0; // end of window (initially the same as start so we can grab the first char)
    int result = 0; // var to store the length of the longest substring we have found so far
    HashSet<Character> chars = new HashSet(); // hashset to track the chars in the current substring
// the right hand side of the window must not go over the bounds of the string
    while(right < s.length()) {
        // use the right hand side of the window to get the char, it is initialised to 0 so will pick the first character first
        if(!chars.contains(s.charAt(right))) {
            chars.add(s.charAt(right));
    // if the current set of unique chars is greater in size than the result, that's our longest substring
            result = Math.max(result, chars.size());
    // increment right to increase the window size and grab a new char
            right++;
        } else {
  // once we encounter a char which is already in the hashset we need create a new window starting from left + 1, and reset right to equal left to initialise a new window
            left ++;
            right = left;
    // reset our tracker
            chars.clear();
        }
    }
    return result;
}
```

### 7. Intersection of 2 Arrays - https://leetcode.com/problems/intersection-of-two-arrays

My solution:
```java
class Solution {
    public int[] intersection(int[] nums1, int[] nums2) {
        HashSet<Integer> nums1Set = new HashSet<Integer>();
        HashSet<Integer> solSet = new HashSet<Integer>();
        
        for(int i = 0; i < nums1.length; i++){
            nums1Set.add(nums1[i]);
        }
        
        for(int i = 0; i < nums2.length; i++){
            if(nums1Set.contains(nums2[i]) && !solSet.contains(nums2[i])){
                solSet.add(nums2[i]);
            }
        }
        
        Iterator<Integer> iter = solSet.iterator();
        int[] solArray = new int[solSet.size()];
        int index = 0;
        
        while(iter.hasNext()){
            solArray[index++] = iter.next();
        }
        
        return solArray;
    }
}
```

### 8. Valid Parentheses - https://leetcode.com/problems/valid-parentheses/

Stacks:
```java
class Solution {
    public boolean isValid(String s) {
        Stack<Character> stack = new Stack<Character>();
        
        for(int i = 0; i < s.length(); i++){
            switch(s.charAt(i)){
                case '(':
                    stack.push('(');
                    break;
                case '{':
                    stack.push('{');
                    break;
                case '[':
                    stack.push('[');
                    break;
                case ')':
                    if(stack.empty() || stack.peek() != '('){
                        return false;
                    } else {
                        stack.pop();
                    }
                    break;
                case '}':
                    if(stack.empty() || stack.peek() != '{'){
                        return false;
                    } else {
                        stack.pop();
                    }
                    break;
                case ']':
                    if(stack.empty() || stack.peek() != '['){
                        return false;
                    } else {
                        stack.pop();
                    }
                    break;
            }
        }
        
        // Are there any lonely brackets which have not been popped?
        if(!stack.empty()){
            return false;
        }
        
        return true;
    }
}
```

### 9. Fibonacci Number - https://leetcode.com/problems/fibonacci-number

Recursion:
```java
class Solution 
{
    public int fib(int N)
    {
        if(N <= 1)
            return N;
        else
            return fib(N - 1) + fib(N - 2);
    }
}
```

### 10. Climbing Stairs - https://leetcode.com/problems/climbing-stairs/

My solution but it exceeds the time limit but still works for most test cases:
```java
class Solution {
   private int result = 0;
  
   public int climbStairs(int n) {
       climb(0, n);
       return this.result;
   }
 
   private void climb(int level, int target) {
       if(level > target) {
           // We have overshot our level, this is useless
           return;
       }
 
       if(level == target){
           this.result++;
           return;
       }
 
       // At every level, we can take 2 decisions - 1 step or 2 steps
       climb(level + 1, target);
       climb(level + 2, target);
   }
}
```

I frankly did not understand this solution very well:
```java
class Solution {
    public int climbStairs(int n) {
        if(n == 0) return 0;
        if(n == 1) return 1;
        if(n == 2) return 2;

        int two_steps_before = 1;
        int one_step_before = 2;
        int all_ways = 0;

        for(int i = 2; i < n; i++){
            all_ways = two_steps_before + one_step_before;
            two_steps_before = one_step_before;
            one_step_before = all_ways;
        }

        return all_ways;
    }
}
```

### 11. Preorder Traversal - https://leetcode.com/problems/binary-tree-preorder-traversal/description/

Trees:
```java
// recursively
public List<Integer> preorderTraversal1(TreeNode root) {
    List<Integer> ret = new ArrayList<>();
    dfs(root, ret);
    return ret;
}

private void dfs(TreeNode root, List<Integer> ret) {
    if (root != null) {
        ret.add(root.val);
        dfs(root.left, ret);
        dfs(root.right, ret);
    }
}
```

### 12. Max Depth of Binary Tree - https://leetcode.com/problems/maximum-depth-of-binary-tree/description/

Trees:
```java
class Solution {
    int maxDepth = 0;

    public int maxDepth(TreeNode root) {
        if(root == null) return maxDepth;
        maxDepth++;
        findMaxDepth(root, 1);

        return maxDepth;
    }

    private void findMaxDepth(TreeNode node, int currentLevel){
        maxDepth = Math.max(maxDepth, currentLevel);
        
        if(node.left != null) findMaxDepth(node.left, currentLevel + 1);
        if(node.right != null) findMaxDepth(node.right, currentLevel + 1);
    }
}
```

### 13. Flood Fill - https://leetcode.com/problems/flood-fill/description/
```java
class Solution {
    public int[][] floodFill(int[][] image, int sr, int sc, int color) {
        changePixels(image, sr, sc, image[sr][sc], color);
        return image;
    }

    private void changePixels(int[][] image, int row, int col, int target, int color) {
        // check out of bounds or incorrect pixel value
        if(
            row < 0 ||
            row > image.length - 1 ||
            col < 0 ||
            col > image[0].length - 1 ||
            image[row][col] != target ||
            image[row][col] == color
        ) {
            return;
        }

        image[row][col] = color;
        
        changePixels(image, row - 1, col, target, color); // up
        changePixels(image, row + 1, col, target, color); // bottom
        changePixels(image, row, col - 1, target, color); // left
        changePixels(image, row, col + 1, target, color); // right
    }
}
```

### 14. Binary Tree Inorder Traversal - https://leetcode.com/problems/binary-tree-inorder-traversal/description/
```java
/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode() {}
 *     TreeNode(int val) { this.val = val; }
 *     TreeNode(int val, TreeNode left, TreeNode right) {
 *         this.val = val;
 *         this.left = left;
 *         this.right = right;
 *     }
 * }
 */
class Solution {
    public List<Integer> inorderTraversal(TreeNode root) {
        List<Integer> solList = new LinkedList<Integer>();
        traverse(root, solList);
        return solList;
    }

    private void traverse (TreeNode node, List<Integer> solList) {
        if(node != null) {
            traverse(node.left, solList);
            solList.add(node.val);
            traverse(node.right, solList);
        }
    }
}
```

### 15. Binary Tree Postorder Traversal - https://leetcode.com/problems/binary-tree-postorder-traversal/description/
```java
/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode() {}
 *     TreeNode(int val) { this.val = val; }
 *     TreeNode(int val, TreeNode left, TreeNode right) {
 *         this.val = val;
 *         this.left = left;
 *         this.right = right;
 *     }
 * }
 */
class Solution {
    public List<Integer> postorderTraversal(TreeNode root) {
        List<Integer> solList = new LinkedList<Integer>();
        traverse(root, solList);
        return solList;
    }

    private void traverse(TreeNode node, List<Integer> solList) {
        if(node != null){
            traverse(node.left, solList);
            traverse(node.right, solList);
            solList.add(node.val);
        }
    }
}
```

### 16. Search in a Binary Search Tree = https://leetcode.com/problems/search-in-a-binary-search-tree/description/
```java
/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode() {}
 *     TreeNode(int val) { this.val = val; }
 *     TreeNode(int val, TreeNode left, TreeNode right) {
 *         this.val = val;
 *         this.left = left;
 *         this.right = right;
 *     }
 * }
 */
class Solution {
    public TreeNode searchBST(TreeNode root, int val) {
        return findNode(root, val);
    }

    private TreeNode findNode(TreeNode node, int val) {
        if(node != null){
            if(node.val == val){
                return node;
            } else if(val > node.val) {
                // Need to search right subtree
                return findNode(node.right, val);
            } else {
                // Need to search left subtree
                return findNode(node.left, val);
            }
        }
        return node;
    }
}
```

### 17. Two Sum - https://leetcode.com/problems/two-sum/description/
```java
import java.util.Map;

class Solution {
    public int[] twoSum(int[] nums, int target) {
        int[] sol = new int[2];
        Map<Integer, Integer> hashMap = new HashMap<>();
        
        for(int i = 0; i < nums.length; i++){
            hashMap.put(nums[i], i);
        }
        
        for(int i = 0; i < nums.length; i++){
            int toBeFound = target - nums[i];
            if(hashMap.containsKey(toBeFound)){
                int solIndex = hashMap.get(toBeFound);
                if(solIndex != i){
                    sol[0] = i;
                    sol[1] = solIndex;
                    break;
                }
            }
        }
        return sol;
    }
}
```

### 18. Remove Duplicates in Sorted Array - https://leetcode.com/problems/remove-duplicates-from-sorted-array/description/

2 pointers:
```java
class Solution {
    public int removeDuplicates(int[] A) {
       if(A.length == 0) return 0;

       int j = 0; // slow pointer

       for(int i = 0; i < A.length; i++) { // fast pointer
           if(A[i] != A[j]) A[++j] = A[i];
       }

       return j + 1;
    }
}
```

### 19. Linked List Cycle - https://leetcode.com/problems/linked-list-cycle/description/
```java
/**
 * Definition for singly-linked list.
 * class ListNode {
 *     int val;
 *     ListNode next;
 *     ListNode(int x) {
 *         val = x;
 *         next = null;
 *     }
 * }
 */
public class Solution {
    public boolean hasCycle(ListNode head) {
        if(head == null) return false;

        HashSet<ListNode> set = new HashSet<>();

        ListNode currentNode = head;
        while(currentNode.next != null){
            set.add(currentNode);
            if(set.contains(currentNode.next)){
                // We have visited the next node before, cycle!
                return true;
            }
            currentNode = currentNode.next;
        }
        return false;
    }
}
```

### Contains Duplicate - https://leetcode.com/problems/contains-duplicate/description/
```java
class Solution {
    public boolean containsDuplicate(int[] nums) {
        HashSet<Integer> set = new HashSet<>();

        for(int i = 0; i < nums.length; i++){
            if(set.contains(nums[i])){
                // We have seen this before, duplicate!
                return true;
            }
            set.add(nums[i]);
        }
        return false;
    }
}
```
