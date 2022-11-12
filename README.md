# Leetcode DSA in Java

## EASY

### Roman to Integer - https://leetcode.com/problems/roman-to-integer
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

### Remove Duplicates in Sorted Array - https://leetcode.com/problems/remove-duplicates-from-sorted-array/description/


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

### Delete Node in Linked List - https://leetcode.com/problems/delete-node-in-a-linked-list

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

### Reverse Linked List - https://leetcode.com/problems/reverse-linked-list/description/

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

### Binary Search - https://leetcode.com/problems/binary-search

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

### Arrange Coins in Stairs - https://leetcode.com/problems/arranging-coins/description/

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

### Longest Substring Without Repeating Characters - https://leetcode.com/problems/longest-substring-without-repeating-characters

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

### Intersection of 2 Arrays - https://leetcode.com/problems/intersection-of-two-arrays

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

### Valid Parentheses - https://leetcode.com/problems/valid-parentheses/

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

### Fibonacci Number - https://leetcode.com/problems/fibonacci-number

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

### Climbing Stairs - https://leetcode.com/problems/climbing-stairs/

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

### Preorder Traversal - https://leetcode.com/problems/binary-tree-preorder-traversal/description/

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

### Max Depth of Binary Tree - https://leetcode.com/problems/maximum-depth-of-binary-tree/description/

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

### Flood Fill - https://leetcode.com/problems/flood-fill/description/
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

### Binary Tree Inorder Traversal - https://leetcode.com/problems/binary-tree-inorder-traversal/description/
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

## MEDIUM

### Delete Nth last node - https://leetcode.com/problems/remove-nth-node-from-end-of-list/description/

You got to remember this only since this has a simple concept but confusing implementation with difficult to handle edge cases.
```java
public ListNode removeNthFromEnd(ListNode head, int n) {
    ListNode start = new ListNode(0);
    ListNode slow = start;
    ListNode fast = start;
    start.next = head;

    for(int i = 0; i <= n; i++ ) {
        fast = fast.next;
    }

    while(fast != null) {
        fast = fast.next;
        slow = slow.next;
    }

    // Delete the node
    slow.next = slow.next.next;

    return start.next;
}
```

### Fruit Into Baskets - https://leetcode.com/problems/fruit-into-baskets/description/ 

My own solution:
```java
class Solution {
    public int totalFruit(int[] fruits) {
        return collectFruitsAroundATree(0, fruits, 0);
    }

    private int collectFruitsAroundATree(int index, int[] fruits, int maxCollected) {
        if(index > fruits.length - 1) return maxCollected;

        int fruitBasket1Type = -1;
        int fruitBasket2Type = -1;
        int currentCollected = 0;

        int i = index;

        // look back first
        while(i >= 0){
            if(fruitBasket1Type == -1 || fruitBasket1Type == fruits[i]){
                fruitBasket1Type = fruits[i--];
                currentCollected++;
            } else if(fruitBasket2Type == -1 || fruitBasket2Type == fruits[i]) {
                fruitBasket2Type = fruits[i--];
                currentCollected++;
            } else {
                // fruit cannot go in any basket!
                break;
            }
        }

        int j = index + 1;

        // look ahead later
        while(j < fruits.length){
            if(fruitBasket1Type == -1 || fruitBasket1Type == fruits[j]){
                fruitBasket1Type = fruits[j++];
                currentCollected++;
            } else if(fruitBasket2Type == -1 || fruitBasket2Type == fruits[j]) {
                fruitBasket2Type = fruits[j++];
                currentCollected++;
            } else {
                // fruit cannot go in any basket!
                break;
            }
        }

        return collectFruitsAroundATree(j, fruits, Math.max(maxCollected, currentCollected));
    }
}
```

### Container with most water - https://leetcode.com/problems/container-with-most-water/description/ 

2 pointers:
```java
class Solution {
    public int maxArea(int[] height) {
        int leftPointer = 0;
        int rightPointer = height.length - 1;
        int maxVolume = 0;

        while(leftPointer < rightPointer) {
            int currentVolume = Math.min(height[leftPointer], height[rightPointer]) * 
            (rightPointer - leftPointer);
            maxVolume = Math.max(maxVolume, currentVolume);
            
            if(height[leftPointer] < height[rightPointer]){
                leftPointer++;
            } else {
                rightPointer--;
            }
        }

        return maxVolume;
    }
}
```

### Implement a Min Stack - https://leetcode.com/problems/min-stack/description/

Stack:
```java
class MinStack {
    LinkedList<Integer> list = new LinkedList<Integer>();
    LinkedList<Integer> minList = new LinkedList<Integer>();
    int size = 0;
    
    public void push(int val) {
        list.push(val);
        size++;
        // we will track the minimum value for each push operation
        minList.push(Math.min(minList.size() == 0 ? Integer.MAX_VALUE : minList.peek(), val));
    }
    
    public void pop() {
        if(size == 0) return;

        list.pop();
        size--;

        minList.pop();
    }
    
    public int top() {
        return list.peek();
    }
    
    public int getMin() {
        return minList.peek();
    }
}

/**
 * Your MinStack object will be instantiated and called as such:
 * MinStack obj = new MinStack();
 * obj.push(val);
 * obj.pop();
 * int param_3 = obj.top();
 * int param_4 = obj.getMin();
 */
```

### Word Search - https://leetcode.com/problems/word-search

Backtracking:
```java
class Solution {
    public boolean exist(char[][] board, String word) {
        for(int i = 0; i < board.length; i++){
            for(int j = 0; j < board[0].length; j++){
                if(explore(board, word, i, j, 0)){
                    return true;   
                }
            }
        }
        
        return false;
    }
    
    private boolean explore(char[][] board, String word, int row, int col, int index){
        if(
            row < 0 ||
            row > board.length - 1 ||
            col < 0 ||
            col > board[0].length - 1 ||
            index > word.length() - 1
        ) {
            // out of bounds
            return false;
        }
        
        if(board[row][col] == word.charAt(index)) {
            // check if this is the last index
            if(index == word.length() - 1){
                return true;
            }
            
            char val = board[row][col];
            board[row][col] = '-';
            
            // Explore all directions since we have a match on this node and
            // return true if any direction matches
            if(
                explore(board, word, row - 1, col, index + 1) || // up
                explore(board, word, row + 1, col, index + 1) || // down
                explore(board, word, row, col - 1, index + 1) || // left
                explore(board, word, row, col + 1, index + 1) // right
            ) {
                return true;
            } else {
                board[row][col] = val;
            }
        }
        
        return false;
    }
}
```

### Rotten Oranges - https://leetcode.com/problems/rotting-oranges/description/
```java
class Solution {
    public int orangesRotting(int[][] grid) {
        int timePassed = -1;
        int freshOranges = 0;
        Queue<int[]> rottenOranges = new LinkedList<int[]>();

        for(int row = 0; row < grid.length; row++){
            for(int col = 0; col < grid[0].length; col++){
                if(grid[row][col] == 2){
                    rottenOranges.add(new int[]{row, col});
                } else if (grid[row][col] == 1) {
                    freshOranges++;
                }
            }
        }

        if(freshOranges == 0) return 0; // All rotten
        if(rottenOranges.size() == 0) return -1; // All fresh

        while(rottenOranges.size() != 0) {
            timePassed++;
            int size = rottenOranges.size();

            for(int i = 0; i < size; i++){
                int[] currentOrange = rottenOranges.remove();
                int row = currentOrange[0];
                int col = currentOrange[1];

                // All surrounding will rot, add to queue
                if(row - 1 >= 0 && grid[row - 1][col] == 1){
                    rottenOranges.add(new int[]{row - 1, col});
                    grid[row - 1][col] = 2; // rotten
                    freshOranges--;
                }
                if(row + 1 < grid.length && grid[row + 1][col] == 1){
                    rottenOranges.add(new int[]{row + 1, col});
                    grid[row + 1][col] = 2; // rotten
                    freshOranges--;
                }
                if(col - 1 >= 0 && grid[row][col - 1] == 1){
                    rottenOranges.add(new int[]{row, col - 1});
                    grid[row][col - 1] = 2; // rotten
                    freshOranges--;
                }
                if(col + 1 < grid[0].length && grid[row][col + 1] == 1){
                    rottenOranges.add(new int[]{row, col + 1});
                    grid[row][col + 1] = 2; // rotten
                    freshOranges--;
                }
            }
        }

        return freshOranges > 0 ? -1 : timePassed;
    }
}
```

### Combination Sum (Can Reuse Element) - https://leetcode.com/problems/combination-sum/description/
```java
class Solution {
    private List<List<Integer>> resultList = new LinkedList<List<Integer>>();
    
    public List<List<Integer>> combinationSum(int[] candidates, int target) {
        Arrays.sort(candidates);
        generateCombinations(0, 0, target, 0, new LinkedList<Integer>(), candidates);

        return resultList;
    }

    private void generateCombinations(
        int numCombs,
        int totalSum,
        int target,
        int start,
        LinkedList<Integer> combs,
        int[] candidates
    ){
        if(numCombs >= 150 || totalSum > target) return;
        if(totalSum == target) this.resultList.add(new LinkedList<Integer>(combs));

        // At any stage, all candidates should be considered, repeat candidates are allowed
        for(int i = start; i < candidates.length; i++) {
            combs.add(candidates[i]);

            generateCombinations(
                numCombs + 1,
                totalSum + candidates[i],
                target,
                i, // not i + 1 since we are also considering duplicates
                combs,
                candidates
            );

            combs.removeLast();
        }
    }
}
```

### Combination Sum (Can’t Reuse Element) - https://leetcode.com/problems/combination-sum-ii/description/
```java
class Solution {
    private List<List<Integer>> resultList = new LinkedList<List<Integer>>();

    public List<List<Integer>> combinationSum2(int[] candidates, int target) {
        Arrays.sort(candidates);
        generateCombs(
            candidates,
            target,
            0,
            0,
            new LinkedList<Integer>()
        );
        return this.resultList;
    }

    private void generateCombs(
        int[] candidates,
        int target,
        int currentSum,
        int start,
        LinkedList<Integer> tempList
    ){
        if(currentSum > target) return; // invalid situation
        if(currentSum == target) this.resultList.add(new LinkedList<Integer>(tempList));

        for(int i = start; i < candidates.length; i++) {
            if(i > start && candidates[i] == candidates[i - 1]) continue; // No duplicates

            tempList.add(candidates[i]);

            generateCombs(
                candidates,
                target,
                currentSum + candidates[i],
                i + 1,
                tempList
            );

            tempList.removeLast();
        }
    }
}
```

### Permutations of Number Array (Without Duplicates) - https://leetcode.com/problems/permutations/description/
```java
class Solution {
    List<List<Integer>> resultList = new LinkedList<List<Integer>>();

    public List<List<Integer>> permute(int[] nums) {
        generatePerms(nums, new LinkedList<Integer>());

        return this.resultList;
    }

    private void generatePerms(int[] nums, LinkedList<Integer> tempList){
        if(tempList.size() == nums.length){
            this.resultList.add(new LinkedList<Integer>(tempList));
            return;
        }

        // At any position we have to consider all values that are left except the values
        // chosen already

        for(int i = 0; i < nums.length; i++) {
            if(tempList.contains(nums[i])) continue; // We have already chosen this element

            tempList.add(nums[i]);
            generatePerms(nums, tempList);
            tempList.removeLast(); // Backtrack cleanup step
        }
    }
}
```

### Permutations of Number Array (With Duplicates) - https://leetcode.com/problems/permutations-ii/description/
```java
class Solution {
    List<List<Integer>> resultList = new LinkedList<List<Integer>>();

    public List<List<Integer>> permuteUnique(int[] nums) {
        Arrays.sort(nums);
        generatePerms(nums, new LinkedList<Integer>(), new boolean[nums.length]);
        return resultList;
    }

    private void generatePerms(
        int[] nums,
        LinkedList<Integer> result,
        boolean[] used
        ) {
        if(result.size() == nums.length){
            this.resultList.add(new LinkedList(result));
            return; // we are done now, backtrack
        }

        for(int i = 0; i < nums.length; i++){
            if(used[i] || (i > 0 && nums[i] == nums[i - 1] && !used[i - 1])) continue;

            result.add(nums[i]);
            used[i] = true;
            generatePerms(nums, result, used);

            // backtrack cleanup step
            result.removeLast();
            used[i] = false;
        }
    }
}
```

### Subsets of Number Array (Without Duplicates) - https://leetcode.com/problems/subsets/description/
```java
class Solution {
    private List<List<Integer>> resultList = new LinkedList<List<Integer>>();

    public List<List<Integer>> subsets(int[] nums) {
        Arrays.sort(nums);
        generateSubsets(nums, new LinkedList<Integer>(), 0);
        return this.resultList;
    }

    private void generateSubsets(int[] nums, LinkedList<Integer> result, int start) {
        if(result.size() > nums.length) {
            return;
        }
        
        this.resultList.add(new LinkedList<Integer>(result));

        for(int i = start; i < nums.length; i++) {
            result.add(nums[i]);
            generateSubsets(nums, result, i + 1);
            result.removeLast();
        }
    }
}
```

### Subsets of Number Array (With Duplicates) - https://leetcode.com/problems/subsets-ii/description/
```java
class Solution {
    List<List<Integer>> resultList = new LinkedList<List<Integer>>();

    public List<List<Integer>> subsetsWithDup(int[] nums) {
        Arrays.sort(nums);
        generateSubsets(nums, new LinkedList<Integer>(), 0);
        return this.resultList;
    }

    private void generateSubsets(
        int[] nums,
        LinkedList<Integer> result,
        int start
    ) {
        this.resultList.add(new LinkedList<Integer>(result));

        for(int i = start; i < nums.length; i++) {
            // This is to prevent duplicates
            if(i > start && nums[i] == nums[i - 1]) continue; 

            result.add(nums[i]);
            generateSubsets(nums, result, i + 1);

            // Backtrack cleanup step
            result.removeLast();
        }
    }
}
```

### Palindrome Partitioning - https://leetcode.com/problems/palindrome-partitioning/description/ 
```java
class Solution {
    List<List<String>> resultList = new LinkedList<List<String>>();

    public List<List<String>> partition(String s) {
        generatePartitions(s, new LinkedList<String>(), 0);

        return this.resultList;
    }

    private void generatePartitions(
        String s,
        LinkedList<String> result,
        int start
    ) {
        if(start == s.length()){
            this.resultList.add(new LinkedList<String>(result));
        } else {
            for(int i = start; i < s.length(); i++) {
                if(this.checkPalindrome(s, start, i)) {
                    result.add(s.substring(start, i + 1));
                    generatePartitions(s, result, i + 1);

                    // backtrack cleanup
                    result.removeLast();
                }
            }
        }
    }

    private boolean checkPalindrome(String s, int low, int high) {
        while(low < high){
            if(s.charAt(low++) != s.charAt(high--)) return false;
        }
        return true;
    }
}
```

### Binary Tree Level Order Traversal - https://leetcode.com/problems/binary-tree-level-order-traversal/description/
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
    List<List<Integer>> solList = new ArrayList<List<Integer>>();

    public List<List<Integer>> levelOrder(TreeNode root) {
        if(root == null) return this.solList;

        traverse(root, 1);
        return this.solList;
    }

    public void traverse(TreeNode node, int level){
        if(this.solList.size() < level){
            // this level has not been reached yet
            this.solList.add(new LinkedList<Integer>());
        }
        this.solList.get(level - 1).add(node.val);

        if(node.left != null) traverse(node.left, level + 1);
        if(node.right != null) traverse(node.right, level + 1);
    }
}
```
