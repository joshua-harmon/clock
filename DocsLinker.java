import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.security.MessageDigest;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * DocsLinker: keeps the "Documents" links at the end of the menu page of tzolkin_clock.html
 * in step with the files in the docs folder next to it.
 *
 * Where the links live in the html: search for "DOCS-START". Between the line  // DOCS-START
 * and the line  // DOCS-END  is a plain list, one line per document:
 *     ['title shown on the menu', 'docs/file name.docx'],
 * You can edit those lines by hand in any text editor, or run this program to rebuild them.
 *
 * What this program does:
 *   1. Reads every file in the docs folder (skipping hidden files, desktop.ini, Thumbs.db and
 *      Word's ~$ lock files).
 *   2. Keeps the title of any file already listed (so titles you changed by hand stay), and makes
 *      a title from the file name for a new file. Files no longer in the folder are dropped.
 *   3. Saves the old html as tzolkin_clock.html.bak, writes the new list between the markers,
 *      changing nothing else, and prints the new size and SHA-256.
 *   4. With -sha, also appends the size and SHA-256 to the newest tzolkin_clock_sha256_*.txt.
 *
 * Editing a document itself (in Word, a text editor, etc.) needs no change to the html: the menu only
 * links to the file. Run this program only after adding, removing or renaming files in docs.
 *
 * Compile and run (Java 8 or later), in the folder that holds tzolkin_clock.html:
 *     javac DocsLinker.java
 *     java DocsLinker                       (uses tzolkin_clock.html and its docs folder)
 *     java DocsLinker -list                 (only shows what is listed now and what is in docs)
 *     java DocsLinker -sha                  (also appends to the SHA-256 record)
 *     java DocsLinker path\to\tzolkin_clock.html
 */
public class DocsLinker {
    static final String START = "// DOCS-START", END = "// DOCS-END";
    static final Pattern ENTRY = Pattern.compile("\\[\\s*'((?:[^'\\\\]|\\\\.)*)'\\s*,\\s*'((?:[^'\\\\]|\\\\.)*)'\\s*\\]");

    public static void main(String[] args) throws Exception {
        boolean listOnly = false, sha = false;
        String htmlPath = "tzolkin_clock.html";
        for (String a : args) {
            if (a.equalsIgnoreCase("-list")) listOnly = true;
            else if (a.equalsIgnoreCase("-sha")) sha = true;
            else htmlPath = a;
        }
        File html = new File(htmlPath).getAbsoluteFile();
        if (!html.isFile()) { System.out.println("Not found: " + html); return; }
        File docsDir = new File(html.getParentFile(), "docs");
        if (!docsDir.isDirectory()) { System.out.println("No docs folder next to the html: " + docsDir); return; }

        String text = new String(Files.readAllBytes(html.toPath()), StandardCharsets.UTF_8);
        String nl = text.contains("\r\n") ? "\r\n" : "\n";
        int s = findMarker(text, START), e = findMarker(text, END);
        if (s < 0 || e < 0 || e < s) { System.out.println("Marker lines " + START + " / " + END + " not found in " + html.getName()); return; }
        int bodyStart = text.indexOf('\n', s) + 1;   // first line after // DOCS-START
        String body = text.substring(bodyStart, e);

        // what is listed now: path -> title, in order
        LinkedHashMap<String, String> listed = new LinkedHashMap<String, String>();
        Matcher m = ENTRY.matcher(body);
        while (m.find()) listed.put(unesc(m.group(2)), unesc(m.group(1)));

        // what is in the docs folder
        File[] fs = docsDir.listFiles();
        List<String> files = new ArrayList<String>();
        if (fs != null) for (File f : fs) {
            String n = f.getName();
            if (!f.isFile() || n.startsWith(".") || n.startsWith("~$") || n.equalsIgnoreCase("desktop.ini") || n.equalsIgnoreCase("Thumbs.db")) continue;
            files.add(n);
        }
        Collections.sort(files, String.CASE_INSENSITIVE_ORDER);

        System.out.println("Listed in " + html.getName() + ":");
        for (Map.Entry<String, String> x : listed.entrySet())
            System.out.println("  " + x.getValue() + "  ->  " + x.getKey() + (new File(html.getParentFile(), x.getKey()).isFile() ? "" : "   (file missing)"));
        System.out.println("In " + docsDir + ":");
        for (String n : files) System.out.println("  " + n + (listed.containsKey("docs/" + n) ? "" : "   (new)"));
        if (listOnly) return;

        StringBuilder nb = new StringBuilder("const DOCS=[" + nl);
        for (String n : files) {
            String p = "docs/" + n, t = listed.containsKey(p) ? listed.get(p) : titleOf(n);
            nb.append(" ['").append(jsEsc(t)).append("','").append(jsEsc(p)).append("'],").append(nl);
        }
        nb.append("];").append(nl);
        String out = text.substring(0, bodyStart) + nb + text.substring(e);
        if (out.equals(text)) { System.out.println("No change: the list already matches the docs folder."); return; }

        File bak = new File(html.getParentFile(), html.getName() + ".bak");
        Files.copy(html.toPath(), bak.toPath(), StandardCopyOption.REPLACE_EXISTING);
        byte[] bytes = out.getBytes(StandardCharsets.UTF_8);
        Files.write(html.toPath(), bytes);
        String hex = sha256(bytes);
        System.out.println("Updated " + html.getName() + " (old copy saved as " + bak.getName() + ")");
        System.out.println("  Size:     " + bytes.length + " bytes");
        System.out.println("  SHA-256:  " + hex);

        if (sha) {
            File rec = newestRecord(html.getParentFile());
            if (rec == null) System.out.println("No tzolkin_clock_sha256_*.txt found; nothing appended.");
            else {
                String line = nl + "New version (documents list updated by DocsLinker, " + new SimpleDateFormat("yyyy-MM-dd HH:mm").format(new Date()) + "):" + nl
                        + "  Size:     " + bytes.length + " bytes" + nl + "  SHA-256:  " + hex + nl;
                Files.write(rec.toPath(), line.getBytes(StandardCharsets.UTF_8), java.nio.file.StandardOpenOption.APPEND);
                System.out.println("Appended to " + rec.getName());
            }
        }
    }

    /** start of the line that is exactly the marker (after trimming), or -1 */
    static int findMarker(String text, String marker) {
        int pos = 0;
        while (pos < text.length()) {
            int eol = text.indexOf('\n', pos);
            if (eol < 0) eol = text.length();
            if (text.substring(pos, eol).trim().equals(marker)) return pos;
            pos = eol + 1;
        }
        return -1;
    }

    /** a title from a file name: drop the extension, underscores to spaces (and hyphens too when there are no spaces) */
    static String titleOf(String name) {
        int d = name.lastIndexOf('.');
        String t = d > 0 ? name.substring(0, d) : name;
        t = t.replace('_', ' ');
        if (!t.contains(" ")) t = t.replace('-', ' ');
        return t.isEmpty() ? name : Character.toUpperCase(t.charAt(0)) + t.substring(1);
    }

    static String jsEsc(String s) { return s.replace("\\", "\\\\").replace("'", "\\'"); }
    static String unesc(String s) { return s.replaceAll("\\\\(.)", "$1"); }

    static File newestRecord(File dir) {
        File best = null;
        File[] fs = dir.listFiles();
        if (fs != null) for (File f : fs)
            if (f.getName().startsWith("tzolkin_clock_sha256_") && f.getName().endsWith(".txt") && (best == null || f.lastModified() > best.lastModified())) best = f;
        return best;
    }

    static String sha256(byte[] b) throws Exception {
        byte[] h = MessageDigest.getInstance("SHA-256").digest(b);
        StringBuilder sb = new StringBuilder();
        for (byte x : h) sb.append(String.format("%02x", x));
        return sb.toString();
    }
}
