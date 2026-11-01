package hypergrass;

import java.util.Random;

/** Round-trip self-test for the order-1 range coder used by auxiliary streams. */
public final class RangeCoderSelfTest {
    public static void main(String[] args) {
        Random rnd = new Random(1234567L);
        check("empty", new byte[0]);
        check("single", new byte[] {42});
        for (int len : new int[] {1, 2, 3, 7, 15, 255, 256, 1000, 65536, 1_000_000}) {
            byte[] random = new byte[len];
            rnd.nextBytes(random);
            check("random-" + len, random);

            byte[] skewed = new byte[len];
            for (int i = 0; i < len; i++) skewed[i] = (byte) (rnd.nextInt(10) == 0 ? 1 : 0);
            check("skewed-" + len, skewed);

            byte[] ascending = new byte[len];
            for (int i = 0; i < len; i++) ascending[i] = (byte) (i & 0xFF);
            check("ascending-" + len, ascending);
        }
        System.out.println("RangeCoder self-test passed");
    }

    private static void check(String name, byte[] data) {
        byte[] comp = RangeCoder.encode(data, data.length);
        byte[] back = RangeCoder.decode(comp, data.length);
        if (back.length != data.length) {
            throw new AssertionError(name + ": length " + back.length + " != " + data.length);
        }
        for (int i = 0; i < data.length; i++) {
            if (back[i] != data[i]) {
                throw new AssertionError(name + ": mismatch at " + i + " ("
                        + (data[i] & 0xFF) + " != " + (back[i] & 0xFF) + ")");
            }
        }
        System.out.println(name + ": " + data.length + " -> " + comp.length);
    }
}
